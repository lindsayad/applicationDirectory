/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    electrostaticFoam

Description
    Solver for electrostatics.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting iteration loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Iteration = " << runTime.timeName() << nl << endl;
		Info<< "Solving equation for the potential." << nl << endl;

        solve
        (
            fvm::laplacian(phi) + eCharge * (pDensity - eDensity) / eps0
        );

		EField = -1.*fvc::grad(phi);
		EMag = mag(EField);
		alpha = alpha0 * exp(-ionCoeff/(EMag + EMagDummy));
		eVelocity = We0 * EField/EMagDummy;
		pVelocity = Wp0 * EField/EMagDummy;
		ePhi = linearInterpolate(eVelocity) & mesh.Sf();
		pPhi = linearInterpolate(pVelocity) & mesh.Sf();
		electronDiffusiveFlux = -De*fvc::grad(eDensity);
		scalarElectronDiffusiveFlux = linearInterpolate(electronDiffusiveFlux) & mesh.Sf()/mag(mesh.Sf());
		scalarIonDiffusiveFlux = linearInterpolate(-Dp*fvc::grad(pDensity)) & mesh.Sf()/mag(mesh.Sf());
		eVelocityDotNormal = linearInterpolate(eVelocity) & mesh.Sf()/mag(mesh.Sf());
		pVelocityDotNormal = linearInterpolate(pVelocity) & mesh.Sf()/mag(mesh.Sf());
		ae = pos(eVelocityDotNormal);
		ap = pos(pVelocityDotNormal);
		pDensityBoundaryCondition = scalarIonDiffusiveFlux/(0.5*vThermalIons + (2.*ap - 2.*unitySurfaceField)*pVelocityDotNormal); 
		gammaIons = pDensity*pVelocity - Dp*fvc::grad(pDensity);
		gammaElectronDensity = (unitySurfaceField - ae) / eMobility * gammaP * mag( (2.*ap - unitySurfaceField) + 0.5 * sqrt( (8. * (massIons + massNeutrals) * massNeutrals) / (pi * (5.*massIons + 3.*massNeutrals) * massIons) ) ) * pMobility * linearInterpolate(pDensity);
		eDensityBoundaryCondition = (-scalarElectronDiffusiveFlux - 0.5*vThermalElectrons*gammaElectronDensity - 2.0*(unitySurfaceField - ae)*gammaP*(linearInterpolate(gammaIons) & mesh.Sf()/mag(mesh.Sf())))/((2.0*unitySurfaceField - 2.0*ae) * eVelocityDotNormal - 0.5*vThermalElectrons);

		Info << "Solving equation for electron density." << nl << endl;

        solve
        (
            fvm::ddt(eDensity)
//		  + fvm::div(ePhi, eDensity)
		  - fvm::laplacian(De, eDensity)
		  - eDensity * alpha * mag(eVelocity)		
        );

		Info << "Solving equation for ion density." << nl << endl;

		solve
		(
            fvm::ddt(pDensity)
//		  + fvm::div(pPhi, pDensity)
		  - fvm::laplacian(Dp, pDensity)
		  - eDensity * alpha * mag(eVelocity)	
		);

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
