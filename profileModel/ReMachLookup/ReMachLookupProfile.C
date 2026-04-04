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

\*---------------------------------------------------------------------------*/

#include "ReMachLookupProfile.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"
#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ReMachLookupProfile, 0);
    addToRunTimeSelectionTable(profileModel, ReMachLookupProfile, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::ReMachLookupProfile::interpolateWeights
(
    const scalar& xIn,
    const List<scalar>& values,
    label& i1,
    label& i2,
    scalar& ddx
) const
{
    // Verbatim copy of lookupProfile::interpolateWeights
    i2 = 0;
    const label nElem = values.size();

    if (nElem == 1)
    {
        i1 = i2;
        ddx = 0.0;
        return;
    }
    else
    {
        while ((i2 < nElem) && (values[i2] < xIn))
        {
            i2++;
        }

        if (i2 == 0)
        {
            i1 = i2;
            ddx = 0.0;
            return;
        }
        else if (i2 == nElem)
        {
            i2 = nElem - 1;
            i1 = i2;
            ddx = 0.0;
            return;
        }
        else
        {
            i1 = i2 - 1;
            ddx = (xIn - values[i1])/(values[i2] - values[i1]);
        }
    }
}


void Foam::ReMachLookupProfile::interpAOA
(
    const label tableIdx,
    const scalar alpha,
    scalar& Cd,
    scalar& Cl
) const
{
    // Mirrors lookupProfile::Cdl exactly — same variable names and pattern
    label i1 = -1;
    label i2 = -1;
    scalar invAlpha = -1.0;
    interpolateWeights(alpha, AOA_[tableIdx], i1, i2, invAlpha);

    Cd = invAlpha*(Cd_[tableIdx][i2] - Cd_[tableIdx][i1]) + Cd_[tableIdx][i1];
    Cl = invAlpha*(Cl_[tableIdx][i2] - Cl_[tableIdx][i1]) + Cl_[tableIdx][i1];
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ReMachLookupProfile::ReMachLookupProfile
(
    const dictionary& dict,
    const word& modelName
)
:
    profileModel(dict, modelName),
    Re_(),
    Mach_(),
    AOA_(),
    Cd_(),
    Cl_()
{
    dict.readEntry("Re",   Re_);
    dict.readEntry("Mach", Mach_);

    const label nRe     = Re_.size();
    const label nMach   = Mach_.size();
    const label nTables = nRe * nMach;

    if (nRe < 1 || nMach < 1)
    {
        FatalIOErrorInFunction(dict)
            << "Re and Mach lists must each contain at least one entry"
            << exit(FatalIOError);
    }

    // Read all tables as a list of (AOA_deg Cd Cl) point lists —
    // same List<vector> element type as lookupProfile uses for its "data"
    List<List<vector>> rawTables;
    dict.readEntry("tables", rawTables);

    if (rawTables.size() != nTables)
    {
        FatalIOErrorInFunction(dict)
            << "Expected " << nTables << " polar tables "
            << "(nRe=" << nRe << " x nMach=" << nMach << ") "
            << "but found " << rawTables.size()
            << exit(FatalIOError);
    }

    AOA_.setSize(nTables);
    Cd_.setSize(nTables);
    Cl_.setSize(nTables);

    forAll(rawTables, ti)
    {
        const List<vector>& raw = rawTables[ti];

        if (raw.empty())
        {
            FatalIOErrorInFunction(dict)
                << "Polar table " << ti << " is empty"
                << exit(FatalIOError);
        }

        AOA_[ti].setSize(raw.size());
        Cd_[ti].setSize(raw.size());
        Cl_[ti].setSize(raw.size());

        // Same loading pattern as lookupProfile constructor
        forAll(raw, i)
        {
            AOA_[ti][i] = degToRad(raw[i][0]);
            Cd_[ti][i]  = raw[i][1];
            Cl_[ti][i]  = raw[i][2];
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ReMachLookupProfile::Cdl
(
    const scalar alpha,
    scalar& Cd,
    scalar& Cl
) const
{
    // Fallback: minimum Re and Mach (table index 0)
    interpAOA(0, alpha, Cd, Cl);
}


void Foam::ReMachLookupProfile::Cdl
(
    const scalar alpha,
    const scalar Re,
    const scalar Ma,
    scalar& Cd,
    scalar& Cl
) const
{
    const label nMach = Mach_.size();

    // Bracket Re
    label iRe1 = -1, iRe2 = -1;
    scalar wRe = 0.0;
    interpolateWeights(Re, Re_, iRe1, iRe2, wRe);

    // Bracket Mach
    label iMa1 = -1, iMa2 = -1;
    scalar wMa = 0.0;
    interpolateWeights(Ma, Mach_, iMa1, iMa2, wMa);

    // Four bounding corner table indices (Re-major storage)
    const label t11 = iRe1*nMach + iMa1;   // (Re1, Ma1)
    const label t12 = iRe1*nMach + iMa2;   // (Re1, Ma2)
    const label t21 = iRe2*nMach + iMa1;   // (Re2, Ma1)
    const label t22 = iRe2*nMach + iMa2;   // (Re2, Ma2)

    // AOA interpolation at each corner — each call mirrors lookupProfile::Cdl
    scalar Cd11, Cl11, Cd12, Cl12, Cd21, Cl21, Cd22, Cl22;
    interpAOA(t11, alpha, Cd11, Cl11);
    interpAOA(t12, alpha, Cd12, Cl12);
    interpAOA(t21, alpha, Cd21, Cl21);
    interpAOA(t22, alpha, Cd22, Cl22);

    // Bilinear combination: first in Mach, then in Re
    const scalar Cd1 = wMa*(Cd12 - Cd11) + Cd11;
    const scalar Cd2 = wMa*(Cd22 - Cd21) + Cd21;
    Cd = wRe*(Cd2 - Cd1) + Cd1;

    const scalar Cl1 = wMa*(Cl12 - Cl11) + Cl11;
    const scalar Cl2 = wMa*(Cl22 - Cl21) + Cl21;
    Cl = wRe*(Cl2 - Cl1) + Cl1;
}


// ************************************************************************* //
