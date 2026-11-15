/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2022 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "rotorDiskSource.H"
#include "addToRunTimeSelectionTable.H"
#include "trimModel.H"
#include "fvMatrices.H"
#include "geometricOneField.H"
#include "syncTools.H"
#include "globalIndex.H"
#include "unitConversion.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(rotorDiskSource, 0);
        addToRunTimeSelectionTable(option, rotorDiskSource, dictionary);
    }
}


const Foam::Enum
<
    Foam::fv::rotorDiskSource::geometryModeType
>
Foam::fv::rotorDiskSource::geometryModeTypeNames_
({
    { geometryModeType::gmAuto, "auto" },
    { geometryModeType::gmSpecified, "specified" },
});


const Foam::Enum
<
    Foam::fv::rotorDiskSource::inletFlowType
>
Foam::fv::rotorDiskSource::inletFlowTypeNames_
({
    { inletFlowType::ifFixed, "fixed" },
    { inletFlowType::ifSurfaceNormal, "surfaceNormal" },
    { inletFlowType::ifLocal, "local" },
});


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fv::rotorDiskSource::checkData()
{
    // Set inflow type
    switch (selectionMode())
    {
        case smCellSet:
        case smCellZone:
        case smAll:
        {
            // Set the profile ID for each blade section
            profiles_.connectBlades(blade_.profileName(), blade_.profileID());
            switch (inletFlow_)
            {
                case ifFixed:
                {
                    coeffs_.readEntry("inletVelocity", inletVelocity_);
                    break;
                }
                case ifSurfaceNormal:
                {
                    scalar UIn(coeffs_.get<scalar>("inletNormalVelocity"));
                    inletVelocity_ = -coordSys_.e3()*UIn;
                    break;
                }
                case ifLocal:
                {
                    break;
                }
                default:
                {
                    FatalErrorInFunction
                        << "Unknown inlet velocity type" << abort(FatalError);
                }
            }

            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Source cannot be used with '"
                << selectionModeTypeNames_[selectionMode()]
                << "' mode.  Please use one of: " << nl
                << selectionModeTypeNames_[smCellSet] << nl
                << selectionModeTypeNames_[smCellZone] << nl
                << selectionModeTypeNames_[smAll]
                << exit(FatalError);
        }
    }
}


void Foam::fv::rotorDiskSource::setFaceArea(vector& axis, const bool correct)
{
    area_ = 0.0;

    static const scalar tol = 0.8;

    const label nInternalFaces = mesh_.nInternalFaces();
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    const vectorField& Sf = mesh_.Sf();
    const scalarField& magSf = mesh_.magSf();

    vector n = Zero;

    // Calculate cell addressing for selected cells
    labelList cellAddr(mesh_.nCells(), -1);
    labelUIndList(cellAddr, cells_) = identity(cells_.size());
    labelList nbrFaceCellAddr(mesh_.nBoundaryFaces(), -1);
    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];

        if (pp.coupled())
        {
            forAll(pp, i)
            {
                label facei = pp.start() + i;
                label nbrFacei = facei - nInternalFaces;
                label own = mesh_.faceOwner()[facei];
                nbrFaceCellAddr[nbrFacei] = cellAddr[own];
            }
        }
    }

    // Correct for parallel running
    syncTools::swapBoundaryFaceList(mesh_, nbrFaceCellAddr);

    // Add internal field contributions
    for (label facei = 0; facei < nInternalFaces; facei++)
    {
        const label own = cellAddr[mesh_.faceOwner()[facei]];
        const label nbr = cellAddr[mesh_.faceNeighbour()[facei]];

        if ((own != -1) && (nbr == -1))
        {
            vector nf = Sf[facei]/magSf[facei];

            if ((nf & axis) > tol)
            {
                area_[own] += magSf[facei];
                n += Sf[facei];
            }
        }
        else if ((own == -1) && (nbr != -1))
        {
            vector nf = Sf[facei]/magSf[facei];

            if ((-nf & axis) > tol)
            {
                area_[nbr] += magSf[facei];
                n -= Sf[facei];
            }
        }
    }


    // Add boundary contributions
    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];
        const vectorField& Sfp = mesh_.Sf().boundaryField()[patchi];
        const scalarField& magSfp = mesh_.magSf().boundaryField()[patchi];

        if (pp.coupled())
        {
            forAll(pp, j)
            {
                const label facei = pp.start() + j;
                const label own = cellAddr[mesh_.faceOwner()[facei]];
                const label nbr = nbrFaceCellAddr[facei - nInternalFaces];
                const vector nf = Sfp[j]/magSfp[j];

                if ((own != -1) && (nbr == -1) && ((nf & axis) > tol))
                {
                    area_[own] += magSfp[j];
                    n += Sfp[j];
                }
            }
        }
        else
        {
            forAll(pp, j)
            {
                const label facei = pp.start() + j;
                const label own = cellAddr[mesh_.faceOwner()[facei]];
                const vector nf = Sfp[j]/magSfp[j];

                if ((own != -1) && ((nf & axis) > tol))
                {
                    area_[own] += magSfp[j];
                    n += Sfp[j];
                }
            }
        }
    }

    if (correct)
    {
        reduce(n, sumOp<vector>());
        axis = n/mag(n);
    }

    if (debug)
    {
        volScalarField area
        (
            IOobject
            (
                name_ + ":area",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimArea, Zero)
        );
        UIndirectList<scalar>(area.primitiveField(), cells_) = area_;

        Info<< type() << ": " << name_ << " writing field " << area.name()
            << endl;

        area.write();
    }
}


void Foam::fv::rotorDiskSource::setUpstreamCells()
{
    // The rotor axis e3 points against the through-flow (thrust acts along
    // +e3, the fluid is accelerated along -e3), so the inflow side of the
    // disk is the +axis side.  setFaceArea() has already collected exactly
    // those faces in area_, including faces on processor and boundary
    // patches, so the upstream cells are the cells with a non-zero area_.

    forAll(cells_, i)
    {
        isUpstreamBoundary_[i] = (area_[i] > ROOTVSMALL);
    }
}


void Foam::fv::rotorDiskSource::buildColumnMap()
{
    // Flood-fill from the upstream boundary cells through axis-aligned faces
    // to assign every zone cell to the upstream cell of its column.
    // Upstream cells are identified by a global zone index so that a column
    // may continue across processor boundaries.
    // Also accumulates the total cell volume per column.

    const vector axis = coordSys_.e3();
    static const scalar tol = 0.8;

    const globalIndex globalZoneCells(cells_.size());
    zoneCellOffset_ = globalZoneCells.toGlobal(0);
    nZoneCellsGlobal_ = globalZoneCells.totalSize();

    const label nInternalFaces = mesh_.nInternalFaces();
    const labelUList& faceOwner = mesh_.faceOwner();
    const labelUList& faceNeighbour = mesh_.faceNeighbour();
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    const vectorField& Sf = mesh_.Sf();
    const scalarField& magSf = mesh_.magSf();

    // Mesh-cell based addressing
    boolList inZone(mesh_.nCells(), false);
    labelList cellSeed(mesh_.nCells(), -1);

    forAll(cells_, i)
    {
        inZone[cells_[i]] = true;

        // Seed: upstream boundary cells map to themselves
        if (isUpstreamBoundary_[i])
        {
            cellSeed[cells_[i]] = zoneCellOffset_ + i;
        }
    }

    bool changed = true;
    while (changed)
    {
        changed = false;

        // Propagate within this processor
        bool changedLocal = true;
        while (changedLocal)
        {
            changedLocal = false;

            for (label facei = 0; facei < nInternalFaces; facei++)
            {
                const label own = faceOwner[facei];
                const label nbr = faceNeighbour[facei];

                if (!inZone[own] || !inZone[nbr]) continue;

                // Only propagate through axially-aligned faces
                if (mag(Sf[facei] & axis) < tol*magSf[facei]) continue;

                if (cellSeed[own] >= 0 && cellSeed[nbr] < 0)
                {
                    cellSeed[nbr] = cellSeed[own];
                    changedLocal = true;
                }
                else if (cellSeed[nbr] >= 0 && cellSeed[own] < 0)
                {
                    cellSeed[own] = cellSeed[nbr];
                    changedLocal = true;
                }
            }

            if (changedLocal)
            {
                changed = true;
            }
        }

        // Propagate across coupled (processor, cyclic) patches
        labelList nbrSeed;
        syncTools::swapBoundaryCellList(mesh_, cellSeed, nbrSeed);

        forAll(pbm, patchi)
        {
            const polyPatch& pp = pbm[patchi];

            if (!pp.coupled()) continue;

            const vectorField& Sfp = mesh_.Sf().boundaryField()[patchi];
            const scalarField& magSfp = mesh_.magSf().boundaryField()[patchi];

            forAll(pp, j)
            {
                const label facei = pp.start() + j;
                const label own = faceOwner[facei];

                if (!inZone[own] || cellSeed[own] >= 0) continue;

                const label seed = nbrSeed[facei - nInternalFaces];

                if (seed < 0) continue;

                if (mag(Sfp[j] & axis) < tol*magSfp[j]) continue;

                cellSeed[own] = seed;
                changed = true;
            }
        }

        reduce(changed, orOp<bool>());
    }

    // Store per zone cell and accumulate the column volumes over all
    // processors, indexed by global zone index of the upstream cell
    scalarList globalColumnVolume(nZoneCellsGlobal_, Zero);
    label nUpstream = 0;
    label nUnresolved = 0;

    forAll(cells_, i)
    {
        const label seed = cellSeed[cells_[i]];

        columnUpstreamIndex_[i] = seed;

        if (seed >= 0)
        {
            globalColumnVolume[seed] += mesh_.V()[cells_[i]];
        }
        else
        {
            ++nUnresolved;
        }

        if (isUpstreamBoundary_[i])
        {
            ++nUpstream;
        }
    }

    Pstream::listCombineReduce(globalColumnVolume, plusEqOp<scalar>());

    forAll(cells_, i)
    {
        columnVolume_[i] =
        (
            isUpstreamBoundary_[i]
          ? globalColumnVolume[zoneCellOffset_ + i]
          : 0.0
        );
    }

    reduce(nUpstream, sumOp<label>());
    reduce(nUnresolved, sumOp<label>());

    Info<< "    Rotor disk columns   = " << nUpstream
        << " (zone cells = " << nZoneCellsGlobal_
        << ", not linked to a column = " << nUnresolved << ")" << endl;

    // Warn if any zone cells were not reachable from an upstream face.
    if (nUnresolved > 0)
    {
        WarningInFunction
            << "buildColumnMap: " << nUnresolved
            << " zone cell(s) could not be linked to an upstream face." << nl
            << "Check that the disk zone has upstream boundary faces "
            << "aligned with the rotor axis." << endl;
    }
}


void Foam::fv::rotorDiskSource::createCoordinateSystem()
{
    // Construct the local rotor coordinate system
    vector origin(Zero);
    vector axis(Zero);
    vector refDir(Zero);

    geometryModeType gm =
        geometryModeTypeNames_.get("geometryMode", coeffs_);

    switch (gm)
    {
        case gmAuto:
        {
            // Determine rotation origin (cell volume weighted)
            scalar sumV = 0.0;
            const scalarField& V = mesh_.V();
            const vectorField& C = mesh_.C();
            forAll(cells_, i)
            {
                const label celli = cells_[i];
                sumV += V[celli];
                origin += V[celli]*C[celli];
            }
            reduce(origin, sumOp<vector>());
            reduce(sumV, sumOp<scalar>());
            origin /= sumV;

            // Determine first radial vector
            vector dx1(Zero);
            scalar magR = -GREAT;
            forAll(cells_, i)
            {
                const label celli = cells_[i];
                vector test = C[celli] - origin;
                if (mag(test) > magR)
                {
                    dx1 = test;
                    magR = mag(test);
                }
            }
            reduce(dx1, maxMagSqrOp<vector>());
            magR = mag(dx1);

            // Determine second radial vector and cross to determine axis
            forAll(cells_, i)
            {
                const label celli = cells_[i];
                vector dx2 = C[celli] - origin;
                if (mag(dx2) > 0.5*magR)
                {
                    axis = dx1 ^ dx2;
                    if (mag(axis) > SMALL)
                    {
                        break;
                    }
                }
            }
            reduce(axis, maxMagSqrOp<vector>());
            axis.normalise();

            // Correct the axis direction using a point above the rotor
            {
                vector pointAbove(coeffs_.get<vector>("pointAbove"));
                vector dir = pointAbove - origin;
                dir.normalise();
                if ((dir & axis) < 0)
                {
                    axis *= -1.0;
                }
            }

            coeffs_.readEntry("refDirection", refDir);

            // Set the face areas and apply correction to calculated axis
            // e.g. if cellZone is more than a single layer in thickness
            setFaceArea(axis, true);

            break;
        }
        case gmSpecified:
        {
            coeffs_.readEntry("origin", origin);
            coeffs_.readEntry("axis", axis);
            coeffs_.readEntry("refDirection", refDir);

            setFaceArea(axis, false);

            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unknown geometryMode " << geometryModeTypeNames_[gm]
                << ". Available geometry modes include "
                << geometryModeTypeNames_
                << exit(FatalError);
        }
    }

    coordSys_ = coordSystem::cylindrical(origin, axis, refDir);

    const scalar sumArea = gSum(area_);
    const scalar diameter = Foam::sqrt(4.0*sumArea/mathematical::pi);
    Info<< "    Rotor gometry:" << nl
        << "    - disk diameter = " << diameter << nl
        << "    - disk area     = " << sumArea << nl
        << "    - origin        = " << coordSys_.origin() << nl
        << "    - r-axis        = " << coordSys_.e1() << nl
        << "    - psi-axis      = " << coordSys_.e2() << nl
        << "    - z-axis        = " << coordSys_.e3() << endl;
}


void Foam::fv::rotorDiskSource::constructGeometry()
{
    const pointUIndList cc(mesh_.C(), cells_);

    // Optional: for later transform(), invTransform()
    /// Rcyl_.reset(coordSys_.R(cc).ptr());

    forAll(cells_, i)
    {
        vector planarPosition = coordSys_.localPosition(cc[i]);

        // Only process cells with non-zero radius (skip hub-centre cells)
        if (planarPosition.x() > ROOTVSMALL)
        {
            x_[i] = planarPosition;

            // Cache max radius
            rMax_ = max(rMax_, x_[i].x());

            // Swept angle relative to rDir axis [radians] in range 0 -> 2*pi
            scalar psi = x_[i].y();

            // Blade flap angle [radians]
            scalar beta =
                flap_.beta0 - flap_.beta1c*cos(psi) - flap_.beta2s*sin(psi);

            // Determine rotation tensor to convert from planar system into the
            // rotor cone system
            scalar c = cos(beta);
            scalar s = sin(beta);
            Rcone_[i] = tensor(c, 0, -s, 0, 1, 0, s, 0, c);
        }
    }
}


Foam::tmp<Foam::vectorField> Foam::fv::rotorDiskSource::inflowVelocity
(
    const volVectorField& U
) const
{
    switch (inletFlow_)
    {
        case ifFixed:
        case ifSurfaceNormal:
        {
            return tmp<vectorField>::New(mesh_.nCells(), inletVelocity_);

            break;
        }
        case ifLocal:
        {
            return U.primitiveField();

            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unknown inlet flow specification" << abort(FatalError);
        }
    }

    return tmp<vectorField>::New(mesh_.nCells(), Zero);
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::fv::rotorDiskSource::rotorDiskSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh

)
:
    fv::cellSetOption(name, modelType, dict, mesh),
    rhoRef_(1.0),
    diameterRef_(1.0), //NEW LINE
    refVelEta_(0.0),//NEW LINE
    omega_(0.0),
    nBlades_(0),
    inletFlow_(ifLocal),
    inletVelocity_(Zero),
    tipEffect_(1.0),
    flap_(),
    x_(cells_.size(), Zero),
    Rcone_(cells_.size(), I),
    area_(cells_.size(), Zero),
    coordSys_(),
    rMax_(0.0),
    trim_(trimModel::New(*this, coeffs_)),
    blade_(coeffs_.subDict("blade")),
    profiles_(coeffs_.subDict("profiles")),
    nu_(1.5e-5),
    speedOfSound_(340.0),
    diskThickness_(0.0),
    diskArea_(0.0),
    diskVolume_(0.0),
    allowMultipleCellThickness_(false),
    solidityFactor_(cells_.size(), Zero),
    isUpstreamBoundary_(cells_.size(), false),
    columnUpstreamIndex_(cells_.size(), -1),
    columnVolume_(cells_.size(), 0.0),
    zoneCellOffset_(0),
    nZoneCellsGlobal_(0)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::rotorDiskSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    volVectorField force
    (
        IOobject
        (
            name_ + ":rotorForce",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedVector(eqn.dimensions()/dimVolume, Zero)
    );

    // Read the reference density for incompressible flow
    coeffs_.readEntry("rhoRef", rhoRef_);

    const vectorField Uin(inflowVelocity(eqn.psi()));
    trim_->correct(Uin, force);
    calculate(geometricOneField(), Uin, trim_->thetag(), force);

    // Add source to rhs of eqn
    eqn -= force;

    if (mesh_.time().writeTime())
    {
        force.write();
    }
}


void Foam::fv::rotorDiskSource::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    volVectorField force
    (
        IOobject
        (
            name_ + ":rotorForce",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedVector(eqn.dimensions()/dimVolume, Zero)
    );

    const vectorField Uin(inflowVelocity(eqn.psi()));
    trim_->correct(rho, Uin, force);
    calculate(rho, Uin, trim_->thetag(), force);

    // Add source to rhs of eqn
    eqn -= force;

    if (mesh_.time().writeTime())
    {
        force.write();
    }
}


bool Foam::fv::rotorDiskSource::read(const dictionary& dict)
{
    if (fv::cellSetOption::read(dict))
    {
        coeffs_.readEntry("fields", fieldNames_);
        fv::option::resetApplied();

        // Read coordinate system/geometry invariant properties
        omega_ = rpmToRads(coeffs_.get<scalar>("rpm"));

        coeffs_.readEntry("nBlades", nBlades_);

        inletFlowTypeNames_.readEntry("inletFlowType", coeffs_, inletFlow_);

        coeffs_.readEntry("tipEffect", tipEffect_);

        // Reference diameter for calculating prandtl tip loss factor
        coeffs_.lookup("diameterRef") >> diameterRef_; //NEW LINE
        
        // Reference velocity for calculating propeller efficiency
        coeffs_.lookup("refVelEta") >> refVelEta_; //NEW LINE

        // Optional: kinematic viscosity and speed of sound for Re/Mach polars
        coeffs_.readIfPresent("nu", nu_);
        coeffs_.readIfPresent("speedOfSound", speedOfSound_);

        // Optional: spread force across multiple axial cell layers
        // Default false = original single-layer behaviour
        coeffs_.readIfPresent("allowMultipleCellThickness", allowMultipleCellThickness_);

        const dictionary& flapCoeffs(coeffs_.subDict("flapCoeffs"));
        flap_.beta0 = degToRad(flapCoeffs.get<scalar>("beta0"));
        flap_.beta1c = degToRad(flapCoeffs.get<scalar>("beta1c"));
        flap_.beta2s = degToRad(flapCoeffs.get<scalar>("beta2s"));


        // Create coordinate system (also calls setFaceArea internally)
        createCoordinateSystem();

        // Identify upstream boundary cells and build column map.
        // Must run after createCoordinateSystem() so coordSys_.e3() is valid.
        if (allowMultipleCellThickness_)
        {
            setUpstreamCells();
            buildColumnMap();
        }

        // Read coordinate system dependent properties
        checkData();

        constructGeometry();
        calculateGeometry();
        calculateSolidityFactor();

        trim_->read(coeffs_);

        if (debug)
        {
            writeField("thetag", trim_->thetag()(), true);
            writeField("faceArea", area_, true);
        }

        return true;
    }

    return false;
}


void Foam::fv::rotorDiskSource::calculateGeometry()
{
    scalar sumArea   = 0.0;
    scalar sumVolume = 0.0;

    forAll(cells_, i)
    {
        sumArea   += area_[i];
        sumVolume += mesh_.V()[cells_[i]];
    }

    diskArea_      = returnReduce(sumArea,   sumOp<scalar>());
    diskVolume_    = returnReduce(sumVolume, sumOp<scalar>());
    diskThickness_ = diskVolume_ / diskArea_;

    Info<< "    Rotor disk thickness = " << diskThickness_ << nl
        << "    Rotor disk area      = " << diskArea_      << nl
        << "    Rotor disk volume    = " << diskVolume_    << endl;
}


void Foam::fv::rotorDiskSource::calculateSolidityFactor()
{
    solidityFactor_ = 0.0;

    forAll(cells_, i)
    {
        const scalar radius = x_[i].x();

        if (radius < SMALL)
        {
            continue;
        }

        scalar twist = 0.0;
        scalar chord = 0.0;
        label  i1    = -1;
        label  i2    = -1;
        scalar invDr = 0.0;
        blade_.interpolate(radius, twist, chord, i1, i2, invDr);

        if (allowMultipleCellThickness_)
        {
            // BEM runs only for the upstream boundary cell of each column.
            // Force density [N/m3] is normalised by the thickness of this
            // column (column volume / inflow face area) so that integrating
            // over the full column volume gives exactly the blade-element
            // force that belongs to the inflow face area of the column.
            // Interior and downstream cells receive their force via Loop 2
            // in calculate() and keep solidityFactor_ = 0.
            if (isUpstreamBoundary_[i])
            {
                const scalar columnThickness = columnVolume_[i]/area_[i];
                solidityFactor_[i] =
                    chord*nBlades_/(mathematical::twoPi*radius*columnThickness);
            }
        }
        else if (area_[i] > ROOTVSMALL)
        {
            // Original single-layer behaviour: boundary-face cells only.
            const scalar cellThickness = mesh_.V()[cells_[i]] / area_[i];
            solidityFactor_[i] =
                chord*nBlades_/(mathematical::twoPi*radius*cellThickness);
        }
        // else: interior cell — solidityFactor_ stays 0
    }
}


// ************************************************************************* //
