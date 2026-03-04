/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "lookupProfile.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(lookupProfile, 0);
    addToRunTimeSelectionTable(profileModel, lookupProfile, dictionary);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::lookupProfile::interpolateWeights
(
    const scalar& xIn,
    const List<scalar>& values,
    label& i1,
    label& i2,
    scalar& ddx
) const
{
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lookupProfile::lookupProfile
(
    const dictionary& dict,
    const word& modelName
)
:
    profileModel(dict, modelName),
    AOA_(),
    Cd_(),
    Cl_(),
    Re_(),
    Ma_(),
    useReMa_(false),
    alphaScale_(1.0),
    ReScale_(1.0),
    MaScale_(1.0)
{
    List<List<scalar>> data;
    if (readFromFile())
    {
        IFstream is(fName_);
        is  >> data;
    }
    else
    {
        dict.readEntry("data", data);
    }

    if (data.empty())
    {
        FatalIOErrorInFunction(dict)
            << "No profile data specified"
            << exit(FatalIOError);
    }

    const label nCols = data[0].size();

    if ((nCols != 3) && (nCols != 5))
    {
        FatalIOErrorInFunction(dict)
            << "Profile data rows must contain 3 values (AOA Cd Cl) or "
            << "5 values (AOA Re Mach Cd Cl)"
            << exit(FatalIOError);
    }

    useReMa_ = (nCols == 5);

    AOA_.setSize(data.size());
    Cd_.setSize(data.size());
    Cl_.setSize(data.size());

    if (useReMa_)
    {
        Re_.setSize(data.size());
        Ma_.setSize(data.size());
    }

    forAll(data, i)
    {
        if (data[i].size() != nCols)
        {
            FatalIOErrorInFunction(dict)
                << "Inconsistent number of columns in profile data at row "
                << i
                << exit(FatalIOError);
        }

        AOA_[i] = degToRad(data[i][0]);

        if (useReMa_)
        {
            Re_[i] = data[i][1];
            Ma_[i] = data[i][2];
            Cd_[i] = data[i][3];
            Cl_[i] = data[i][4];
        }
        else
        {
            Cd_[i] = data[i][1];
            Cl_[i] = data[i][2];
        }
    }

    if (useReMa_)
    {
        alphaScale_ = max(max(AOA_) - min(AOA_), SMALL);
        ReScale_ = max(max(Re_) - min(Re_), SMALL);
        MaScale_ = max(max(Ma_) - min(Ma_), SMALL);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::lookupProfile::Cdl(const scalar alpha, scalar& Cd, scalar& Cl) const
{
    if (useReMa_)
    {
        Cdl(alpha, 0.0, 0.0, Cd, Cl);
        return;
    }

    label i1 = -1;
    label i2 = -1;
    scalar invAlpha = -1.0;
    interpolateWeights(alpha, AOA_, i1, i2, invAlpha);

    Cd = invAlpha*(Cd_[i2] - Cd_[i1]) + Cd_[i1];
    Cl = invAlpha*(Cl_[i2] - Cl_[i1]) + Cl_[i1];
}


void Foam::lookupProfile::Cdl
(
    const scalar alpha,
    const scalar Re,
    const scalar Ma,
    scalar& Cd,
    scalar& Cl
) const
{
    if (!useReMa_)
    {
        Cdl(alpha, Cd, Cl);
        return;
    }

    scalar sumW = 0.0;
    scalar CdW = 0.0;
    scalar ClW = 0.0;

    forAll(AOA_, i)
    {
        const scalar da = (alpha - AOA_[i])/alphaScale_;
        const scalar dr = (Re - Re_[i])/ReScale_;
        const scalar dm = (Ma - Ma_[i])/MaScale_;

        const scalar d2 = da*da + dr*dr + dm*dm;
        const scalar w = 1.0/max(d2, SMALL);

        sumW += w;
        CdW += w*Cd_[i];
        ClW += w*Cl_[i];
    }

    Cd = CdW/sumW;
    Cl = ClW/sumW;
}


// ************************************************************************* //
