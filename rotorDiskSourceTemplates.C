/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
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
#include "volFields.H"
#include "unitConversion.H"

using namespace Foam::constant;

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class RhoFieldType>
void Foam::fv::rotorDiskSource::calculate
(
    const RhoFieldType& rho,
    const vectorField& U,
    const scalarField& thetag,
    vectorField& force,
    const bool divideVolume,
    const bool output
) const
{
    const scalarField& V = mesh_.V();

    // Logging info
    scalar dragEff = 0.0;
    scalar liftEff = 0.0;
    scalar AOAmin = GREAT;
    scalar AOAmax = -GREAT;

    // Cached position-dependent rotations available?
    const bool hasCache = bool(Rcyl_);

    forAll(cells_, i)
    {
        if (area_[i] > ROOTVSMALL)
        {
            const label celli = cells_[i];

            const scalar radius = x_[i].x();

            const tensor Rcyl =
            (
                hasCache
              ? (*Rcyl_)[i]
              : coordSys_.R(mesh_.C()[celli])
            );

            // Transform velocity into local cylindrical reference frame
            vector Uc = invTransform(Rcyl, U[celli]);

            // Transform velocity into local coning system
            Uc = transform(Rcone_[i], Uc);

            // Set radial component of velocity to zero
            Uc.x() = 0.0;

            // Set blade normal component of velocity
            Uc.y() = radius*omega_ - Uc.y();

            // Determine blade data for this radius
            // i2 = index of upper radius bound data point in blade list
            scalar twist = 0.0;
            scalar chord = 0.0;
            label i1 = -1;
            label i2 = -1;
            scalar invDr = 0.0;
            blade_.interpolate(radius, twist, chord, i1, i2, invDr);

            // Flip geometric angle if blade is spinning in reverse (clockwise)
            scalar alphaGeom = thetag[i] + twist;
            if (omega_ < 0)
            {
                alphaGeom = mathematical::pi - alphaGeom;
            }

            // Effective angle of attack
            scalar alphaEff = alphaGeom - atan2(-Uc.z(), Uc.y());
            if (alphaEff > mathematical::pi)
            {
                alphaEff -= mathematical::twoPi;
            }
            if (alphaEff < -mathematical::pi)
            {
                alphaEff += mathematical::twoPi;
            }

            AOAmin = min(AOAmin, alphaEff);
            AOAmax = max(AOAmax, alphaEff);

            // Determine section profile data for this angle of attack,
            // local Reynolds number and local Mach number without radial
            // interpolation of profile polars.
            const label profilei = blade_.profileID()[invDr < 0.5 ? i1 : i2];

            const scalar Urel = mag(Uc);
            const scalar Re = Urel*chord/max(nu_, SMALL);
            const scalar Ma = Urel/max(aRef_, SMALL);

            scalar Cd = 0.0;
            scalar Cl = 0.0;
            profiles_[profilei].Cdl(alphaEff, Re, Ma, Cd, Cl);

            // Apply tip effect for blade lift
            scalar tipFactor = neg(radius/rMax_ - tipEffect_);

            // Calculate forces perpendicular to blade
            scalar pDyn = 0.5*rho[celli]*magSqr(Uc);

            scalar f = pDyn*chord*nBlades_*area_[i]/radius/mathematical::twoPi;
            vector localForce = vector(0.0, -f*Cd, tipFactor*f*Cl);

            // Accumulate forces
            dragEff += rhoRef_*localForce.y();
            liftEff += rhoRef_*localForce.z();

            // Transform force from local coning system into rotor cylindrical
            localForce = invTransform(Rcone_[i], localForce);

            // Transform force into global Cartesian coordinate system
            force[celli] = transform(Rcyl, localForce);

            if (divideVolume)
            {
                force[celli] /= V[celli];
            }
        }
    }

    if (output)
    {
        reduce(AOAmin, minOp<scalar>());
        reduce(AOAmax, maxOp<scalar>());
        reduce(dragEff, sumOp<scalar>());
        reduce(liftEff, sumOp<scalar>());

        Info<< type() << " output:" << nl
            << "    min/max(AOA)   = " << radToDeg(AOAmin) << ", "
            << radToDeg(AOAmax) << nl
            << "    Effective drag = " << dragEff << nl
            << "    Effective lift = " << liftEff << endl;
    }
}


template<class Type>
void Foam::fv::rotorDiskSource::writeField
(
    const word& name,
    const List<Type>& values,
    const bool writeNow
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> FieldType;

    if (mesh_.time().writeTime() || writeNow)
    {
        auto tfield = tmp<FieldType>::New
        (
            IOobject
            (
                name,
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensioned<Type>(dimless, Zero)
        );

        auto& field = tfield.ref().primitiveFieldRef();

        if (cells_.size() != values.size())
        {
            FatalErrorInFunction
                << abort(FatalError);
        }

        forAll(cells_, i)
        {
            const label celli = cells_[i];
            field[celli] = values[i];
        }

        tfield().write();
    }
}


// ************************************************************************* //
