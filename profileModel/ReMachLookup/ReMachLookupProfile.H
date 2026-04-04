/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
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

Class
    Foam::ReMachLookupProfile

Description
    Profile model where lift and drag coefficients are interpolated from a
    three-dimensional table keyed on Reynolds number, Mach number and angle
    of attack.

    Per-table AOA interpolation uses the same interpolateWeights function
    as lookupProfile (copied verbatim).  Bilinear interpolation is then
    applied across (Re, Mach).

Usage
    \verbatim
    profiles
    {
        myProfile
        {
            type    ReMachLookup;

            // Sorted Reynolds-number breakpoints
            Re      (100000  300000  1000000);

            // Sorted Mach-number breakpoints
            Mach    (0.1  0.2  0.3);

            // Polar tables in Re-major order:
            //   entry[iRe*nMach + iMach] covers (Re[iRe], Mach[iMach])
            // Each row: (AOA_deg  Cd  Cl)
            tables
            (
                // Re=1e5  Ma=0.1
                ((-180 0.020 0.000)(-90 1.800 0.000)(0 0.011 0.000)(90 1.800 0.000)(180 0.020 0.000))
                // Re=1e5  Ma=0.2
                (...)
                ...
            );
        }
    }
    \endverbatim

Note
  - AOA values are supplied in degrees and converted to radians internally.
  - Re and Mach values must be supplied in ascending order.
  - Extrapolation is clamped to the table boundaries.

See also
  - Foam::profileModel
  - Foam::lookupProfile

SourceFiles
    ReMachLookupProfile.C

\*---------------------------------------------------------------------------*/

#ifndef ReMachLookupProfile_H
#define ReMachLookupProfile_H

#include "profileModel.H"
#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class ReMachLookupProfile Declaration
\*---------------------------------------------------------------------------*/

class ReMachLookupProfile
:
    public profileModel
{
    // Private Data

        //- Sorted Reynolds-number breakpoints
        List<scalar> Re_;

        //- Sorted Mach-number breakpoints
        List<scalar> Mach_;

        //- Per-table AOA values [rad], indexed [iRe*nMach + iMach]
        List<List<scalar>> AOA_;

        //- Per-table drag coefficient values
        List<List<scalar>> Cd_;

        //- Per-table lift coefficient values
        List<List<scalar>> Cl_;


    // Private Member Functions

        //- Interpolation weights — verbatim copy of lookupProfile::interpolateWeights
        void interpolateWeights
        (
            const scalar& xIn,
            const List<scalar>& values,
            label& i1,
            label& i2,
            scalar& ddx
        ) const;

        //- Interpolate Cd and Cl in AOA for one (Re,Mach) table.
        //  Mirrors lookupProfile::Cdl exactly.
        void interpAOA
        (
            const label tableIdx,
            const scalar alpha,
            scalar& Cd,
            scalar& Cl
        ) const;


public:

    //- Runtime type information
    TypeName("ReMachLookup");


    // Constructors

        //- Construct from dictionary and model name
        ReMachLookupProfile(const dictionary& dict, const word& modelName);

        //- No copy construct
        ReMachLookupProfile(const ReMachLookupProfile&) = delete;

        //- No copy assignment
        void operator=(const ReMachLookupProfile&) = delete;


    //- Destructor
    ~ReMachLookupProfile() = default;


    // Member Functions

        //- Return Cd and Cl for given AOA only (uses minimum-Re/Mach table)
        virtual void Cdl
        (
            const scalar alpha,
            scalar& Cd,
            scalar& Cl
        ) const;

        //- Return Cd and Cl for given AOA, Reynolds number and Mach number
        virtual void Cdl
        (
            const scalar alpha,
            const scalar Re,
            const scalar Ma,
            scalar& Cd,
            scalar& Cl
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
