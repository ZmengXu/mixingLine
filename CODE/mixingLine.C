/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

Application
    mixingLine

Description
    Get a mixing line for adiabatic mixing of fuel and oxidizer,
    designed for use on single cell cases to calculate the mixed thermo,
    similar to chemFoam

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "rhoReactionThermo.H"
#include "reactingMixture.H"
#include "thermoPhysicsTypes.H"
#include "basicSpecieMixture.H"
#include "thermoTypeFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();

    #define CREATE_MESH createSingleCellMesh.H
    #define NO_CONTROL

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createSingleCellMesh.H"
    #include "createFields.H"
    #include "readInitialConditions.H"

    
    forAll(Z, iRatio)
    {
        const scalar targetMixingRatio = Z[iRatio];
        
        if( targetMixingRatio<0.0 || targetMixingRatio>1.0 )
        {
            FatalError
                << "in initialConditions, targetMixingRatio "
                << targetMixingRatio << nl
                << "    is out of valid range: [0, 1]."
                << exit(FatalError);
        }
        
        const scalar hmix = targetMixingRatio*hfuel+(1.0-targetMixingRatio)*hoxidizer;

        scalarList Ymix(nSpecie, 0.0);
        Ymix = targetMixingRatio*Yfuel+(1.0-targetMixingRatio)*Yoxidizer;
        const scalar mTotMix = sum(Yoxidizer);
        forAll(Y, i)
        {
            Ymix[i] /= mTotMix;
        }
        

        thermo.he() = dimensionedScalar(dimEnergy/dimMass, hmix);
        forAll(Y, i)
        {
            Y[i] = Ymix[i];
        }

        thermo.correct();

        forAll(Y, speciesI)
        {
            writeData[speciesI].append(Y[speciesI][0]);
        }
        writeData[nSpecie].append(thermo.T()[0]);
        writeData[nSpecie+1].append(thermo.rho()[0]);
        
        Info<< " p   = " << p[0] << " [Pa]" << nl
            << " T   = " << thermo.T()[0] << " [K] " << nl
            << " rho   = " << thermo.rho()[0] << " [kg/m3] " << nl
            << " Ymix = " << Ymix << nl
            << endl;
    }

    forAll(writeData, i)
    {
        writeData[i].write();
    }
    
    Info << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
