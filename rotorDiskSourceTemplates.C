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
#include "ReMachLookupProfile.H"

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
    scalar thrustEff = 0.0; //NEW LINE
    scalar torqueEff = 0.0; //NEW LINE
    scalar powerEff = 0.0; //NEW LINE
    scalar AOAmin = GREAT;
    scalar AOAmax = -GREAT;
    scalar epsMin = GREAT; //NEW LINE
    scalar epsMax = -GREAT; //NEW LINE

    // Cached position-dependent rotations available?
    const bool hasCache = bool(Rcyl_);

    forAll(cells_, i)
    {
        if (solidityFactor_[i] < ROOTVSMALL)
        {
            continue;
        }

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

        // Determine profile data for this radius and angle of attack
        const label profile1 = blade_.profileID()[i1];
        const label profile2 = blade_.profileID()[i2];

        // Local flow speed and derived dimensionless numbers
        const scalar Umag = mag(Uc);
        const scalar Re   = (nu_ > VSMALL) ? Umag*chord/nu_ : 0.0;
        const scalar Ma   = (speedOfSound_ > VSMALL) ? Umag/speedOfSound_ : 0.0;

        if (debug)
        {
            Info<< "rotorDiskSource DEBUG cell=" << celli
                << " r=" << radius
                << " alphaEff=" << radToDeg(alphaEff) << " deg"
                << " Umag=" << Umag
                << " Re=" << Re
                << " Ma=" << Ma
                << " chord=" << chord
                << nl;
        }

        scalar Cd1 = 0.0;
        scalar Cl1 = 0.0;
        {
            const ReMachLookupProfile* p =
                dynamic_cast<const ReMachLookupProfile*>(&profiles_[profile1]);
            if (p)
                p->Cdl(alphaEff, Re, Ma, Cd1, Cl1);
            else
                profiles_[profile1].Cdl(alphaEff, Cd1, Cl1);
        }

        scalar Cd2 = 0.0;
        scalar Cl2 = 0.0;
        {
            const ReMachLookupProfile* p =
                dynamic_cast<const ReMachLookupProfile*>(&profiles_[profile2]);
            if (p)
                p->Cdl(alphaEff, Re, Ma, Cd2, Cl2);
            else
                profiles_[profile2].Cdl(alphaEff, Cd2, Cl2);
        }

        scalar Cd = invDr*(Cd2 - Cd1) + Cd1;
        scalar Cl = invDr*(Cl2 - Cl1) + Cl1;

        // Apply tip effect for blade lift
        scalar tipFactor = neg(radius/rMax_ - tipEffect_);

        // Flow angle
        scalar eps = atan2(-Uc.z(), Uc.y());
        if (eps < -mathematical::pi)
        {
            eps = (2.0*mathematical::pi + eps);
        }
        if (eps > mathematical::pi)
        {
            eps = (eps - 2.0*mathematical::pi);
        }
        epsMin = min(epsMin, eps);
        epsMax = max(epsMax, eps);

        //Drela tip factor:
        //scalar lambdaVal = radius/(0.5*diameterRef_)*tan(eps);
        //scalar tipFactor_f = (nBlades_/2)*(1-radius/(0.5*diameterRef_))*(1/lambdaVal);
        //scalar tipFactor = 2/mathematical::pi*acos(exp(-tipFactor_f));

        // Dynamic pressure
        scalar pDynamic = 0.5*rho[celli]*magSqr(Uc);

        // Corrected for distribution and 3D effects
        pDynamic *= solidityFactor_[i]*tipFactor;

        // Tangential and axial forces
        scalar fTangen = pDynamic*( Cd*cos(eps) + Cl*sin(eps));
        scalar fAxial  = pDynamic*(-Cd*sin(eps) + Cl*cos(eps));

        vector localForce = vector(0.0, -fTangen, fAxial);

        // Accumulate forces.
        // In multi-layer mode (allowMultipleCellThickness_) only upstream cells
        // enter this loop, so use the full column volume for integration.
        const scalar integrationVolume =
            allowMultipleCellThickness_ ? columnVolume_[i] : V[celli];

        dragEff += rhoRef_*localForce.y();
        liftEff += rhoRef_*localForce.z();

        thrustEff += rhoRef_*fAxial  *integrationVolume;
        torqueEff += rhoRef_*fTangen *radius*integrationVolume;
        powerEff  += rhoRef_*fTangen *radius*omega_*integrationVolume;

        // Transform force from local coning system into rotor cylindrical
        localForce = invTransform(Rcone_[i], localForce);

        // Transform force into global Cartesian coordinate system
        // pDynamic*solidityFactor_ is already force density [N/m3] — no V division needed
        force[celli] = transform(Rcyl, localForce);
    }

    // Loop 2 - distribute upstream force to interior and downstream cells.
    // Each cell in a column receives the same force density as its upstream
    // boundary cell, which may live on another processor.
    // Only active in multi-layer mode.
    if (allowMultipleCellThickness_)
    {
        List<vector> upstreamForce(nZoneCellsGlobal_, Zero);

        forAll(cells_, i)
        {
            if (isUpstreamBoundary_[i])
            {
                upstreamForce[zoneCellOffset_ + i] = force[cells_[i]];
            }
        }

        Pstream::listCombineReduce(upstreamForce, plusEqOp<vector>());

        forAll(cells_, i)
        {
            const label upIdx = columnUpstreamIndex_[i];
            if (upIdx >= 0 && upIdx != zoneCellOffset_ + i)
            {
                force[cells_[i]] = upstreamForce[upIdx];
            }
        }
    }

    if (output)
    {
        reduce(AOAmin,    minOp<scalar>());
        reduce(AOAmax,    maxOp<scalar>());
        reduce(epsMin,    minOp<scalar>());
        reduce(epsMax,    maxOp<scalar>());
        reduce(dragEff,   sumOp<scalar>());
        reduce(liftEff,   sumOp<scalar>());
        reduce(thrustEff, sumOp<scalar>());
        reduce(torqueEff, sumOp<scalar>());
        reduce(powerEff,  sumOp<scalar>());

        // Efficiency from the globally summed thrust and power
        const scalar etaProp =
            (mag(powerEff) > VSMALL) ? thrustEff*refVelEta_/powerEff : 0.0;

        Info<< type() << " output:" << nl
            << "    min/max(AOA)   = " << radToDeg(AOAmin) << ", "
            << radToDeg(AOAmax) << nl
            << " min/max(eps) = " << radToDeg(epsMin) << ", " //NEW LINE
            << radToDeg(epsMax) << nl //NEW LINE
            << " Rotor thrust = " << thrustEff << nl //NEW LINE
            << " Rotor torque = " << torqueEff << nl //NEW LINE
            << " Rotor power = " << powerEff << nl //NEW LINE
            << " Rotor propeller efficiency = " << etaProp << nl //NEW LINE
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
