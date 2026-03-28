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
#include "vector.H"
#include "unitConversion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ReMachLookupProfile, 0);
    addToRunTimeSelectionTable(profileModel, ReMachLookupProfile, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::ReMachLookupProfile::interpolateWeights
(
    const scalar xIn,
    const List<scalar>& values,
    label& i1,
    label& i2,
    scalar& ddx
) const
{
    const label nElem = values.size();

    if (nElem == 1)
    {
        i1 = 0;
        i2 = 0;
        ddx = 0.0;
        return;
    }

    i2 = 0;
    while (i2 < nElem && values[i2] < xIn)
    {
        i2++;
    }

    if (i2 == 0)
    {
        i1 = 0;
        ddx = 0.0;
        return;
    }

    if (i2 == nElem)
    {
        i1 = nElem - 1;
        i2 = nElem - 1;
        ddx = 0.0;
        return;
    }

    i1 = i2 - 1;
    ddx = (xIn - values[i1])/(values[i2] - values[i1]);
}


void Foam::ReMachLookupProfile::interpAOA
(
    const label tableIdx,
    const scalar alpha,
    scalar& Cd,
    scalar& Cl
) const
{
    const List<scalar>& AOA = AOA_[tableIdx];
    const List<scalar>& CdT = Cd_[tableIdx];
    const List<scalar>& ClT = Cl_[tableIdx];

    label i1 = -1;
    label i2 = -1;
    scalar w = 0.0;
    interpolateWeights(alpha, AOA, i1, i2, w);

    Cd = w*(CdT[i2] - CdT[i1]) + CdT[i1];
    Cl = w*(ClT[i2] - ClT[i1]) + ClT[i1];
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
    dict.readEntry("Re", Re_);
    dict.readEntry("Mach", Mach_);

    const label nRe   = Re_.size();
    const label nMach = Mach_.size();
    const label nTables = nRe * nMach;

    if (nRe < 1 || nMach < 1)
    {
        FatalIOErrorInFunction(dict)
            << "Re and Mach lists must each contain at least one entry"
            << exit(FatalIOError);
    }

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
    // Fallback: use first table entry (minimum Re and Mach)
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

    // Interpolation weights in Re
    label iRe1 = -1, iRe2 = -1;
    scalar wRe = 0.0;
    interpolateWeights(Re, Re_, iRe1, iRe2, wRe);

    // Interpolation weights in Mach
    label iMa1 = -1, iMa2 = -1;
    scalar wMa = 0.0;
    interpolateWeights(Ma, Mach_, iMa1, iMa2, wMa);

    // Table indices for the four bounding corners
    const label t11 = iRe1*nMach + iMa1;
    const label t12 = iRe1*nMach + iMa2;
    const label t21 = iRe2*nMach + iMa1;
    const label t22 = iRe2*nMach + iMa2;

    // AOA interpolation at each corner
    scalar Cd11, Cl11, Cd12, Cl12, Cd21, Cl21, Cd22, Cl22;
    interpAOA(t11, alpha, Cd11, Cl11);
    interpAOA(t12, alpha, Cd12, Cl12);
    interpAOA(t21, alpha, Cd21, Cl21);
    interpAOA(t22, alpha, Cd22, Cl22);

    // Bilinear combination in (Re, Mach)
    const scalar Cd1 = wMa*(Cd12 - Cd11) + Cd11;
    const scalar Cd2 = wMa*(Cd22 - Cd21) + Cd21;
    Cd = wRe*(Cd2 - Cd1) + Cd1;

    const scalar Cl1 = wMa*(Cl12 - Cl11) + Cl11;
    const scalar Cl2 = wMa*(Cl22 - Cl21) + Cl21;
    Cl = wRe*(Cl2 - Cl1) + Cl1;
}


// ************************************************************************* //
