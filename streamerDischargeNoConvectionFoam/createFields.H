    Info<< "Reading physicalProperties\n" << endl;

    IOdictionary physicalProperties
    (
        IOobject
        (
            "physicalProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar alpha0
    (
        physicalProperties.lookup("alpha0")
    );

	dimensionedScalar ionCoeff
	(
		physicalProperties.lookup("ionCoeff")
	);

	dimensionedScalar We0
	(
		physicalProperties.lookup("We0")
	);

	dimensionedScalar WeExp
	(
		physicalProperties.lookup("WeExp")
	);

	dimensionedScalar Wp0
	(
		physicalProperties.lookup("Wp0")
	);

	dimensionedScalar De
	(
		physicalProperties.lookup("De")
	);

	dimensionedScalar Dp
	(
		physicalProperties.lookup("Dp")
	);

	dimensionedScalar eCharge
	(
		physicalProperties.lookup("eCharge")
	);

	dimensionedScalar eps0
	(
		physicalProperties.lookup("eps0")
	);

	dimensionedScalar eMobility
	(
		physicalProperties.lookup("eMobility")
	);
		
	dimensionedScalar pMobility
	(
		physicalProperties.lookup("pMobility")
	);

    Info<< "Reading field phi\n" << endl;
    volScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

	Info<< "Here is the boundary field for phi: " << phi.boundaryField() << nl << endl;
	Info<< "Here is the internal field for phi: " << phi.internalField() << nl << endl;

    Info<< "Reading field eDensity\n" << endl;
    volScalarField eDensity
    (
        IOobject
        (
            "eDensity",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

	Info<< "Reading field pDensity\n" << endl;
	volScalarField pDensity
	(
		IOobject
		(
			"pDensity",
			runTime.timeName(),
			mesh,
			IOobject::MUST_READ,
			IOobject::AUTO_WRITE
		),
		mesh
	);

	Info<< "Calculating field EField\n" << endl;
	volVectorField EField
	(
		IOobject
		(
			"EField",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		-1.*fvc::grad(phi)
	);		

	Info<< "Calculating field EMag\n" << endl;
	volScalarField EMag
	(
		IOobject
		(
			"EMag",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mag(EField)
	);

	Info<< "Defining EMagDummy\n" << endl;
	volScalarField EMagDummy
	(
		IOobject
		(
			"EMagDummy",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mesh,
		dimensionedScalar("EMagDummy", dimensionSet(1, 1, -3, 0, 0, -1, 0), Foam::scalar(1.))
	);
		
	Info<< "Calculating field alpha\n" << endl;
	volScalarField alpha
	(
		IOobject
		(
			"alpha",
			runTime.timeName(),		
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		alpha0*exp(-ionCoeff/(EMag+EMagDummy))
	);

	Info<< "Calculating eVelocity\n" << endl;
	volVectorField eVelocity
	(
		IOobject
		(
			"eVelocity",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		We0 * EField/EMagDummy
//		We0 * (scale(EField,scale(EField,EField))*EMagDummy / scale(EField,scale(EField,scale(EField,EField))))
	);

	Info<< "Calculating pVelocity\n" << endl;
	volVectorField pVelocity
	(
		IOobject
		(
			"pVelocity",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		Wp0 * EField/EMagDummy
	);

	Info<< "Calculating face flux field ePhi\n" << endl;
	surfaceScalarField ePhi
	(
		IOobject
		(
			"ePhi",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		linearInterpolate(eVelocity) & mesh.Sf()
	);

	Info<< "Calculating face flux field pPhi\n" << endl;
	surfaceScalarField pPhi
	(
		IOobject
		(
			"pPhi",
			runTime.timeName(),
			mesh,
			IOobject::READ_IF_PRESENT,
			IOobject::AUTO_WRITE
		),
		linearInterpolate(pVelocity) & mesh.Sf()
	);

/*	Info << "Calculating ae\n" << endl;
	surfaceScalarField ae
	(
		IOobject
		(
			"ae",	
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		((linearInterpolate(eVelocity) & mesh.Sf()) > 0) ? 1. : 0.
	);

	Info << "Calculating ap\n" << endl;
	surfaceScalarField ap
	(
		IOobject
		(
			"ap",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		((linearInterpolate(pVelocity) & mesh.Sf()) > 0) ? 1. : 0.
	); */

	volVectorField electronDiffusiveFlux
	(
		IOobject
		(
			"electronDiffusiveFlux",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		-De*fvc::grad(eDensity)
	);

	surfaceScalarField scalarElectronDiffusiveFlux
	(
		IOobject
		(
			"scalarElectronDiffusiveFlux",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		linearInterpolate(electronDiffusiveFlux) & mesh.Sf()/mag(mesh.Sf())
	);

	surfaceScalarField scalarIonDiffusiveFlux
	(
		IOobject
		(
			"scalarIonDiffusiveFlux",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		linearInterpolate(-Dp*fvc::grad(pDensity)) & mesh.Sf()/mag(mesh.Sf())
	);

	surfaceScalarField scalarElectronConvectiveFluxDensity
	(
		IOobject
		(
			"scalarElectronConvectiveFlux",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		linearInterpolate(eVelocity*eDensity) & mesh.Sf()/mag(mesh.Sf())
	);

	surfaceScalarField eVelocityDotNormal
	(
		IOobject
		(
			"eVelocityDotNormal",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		linearInterpolate(eVelocity) & mesh.Sf()/mag(mesh.Sf())
	);	

	Info<< "Here is eVelocityDotNormal " << nl << eVelocityDotNormal << nl << endl;
	Info<< "Here is pos(eVelocityDotNormal) " << nl << pos(eVelocityDotNormal) << nl << endl;

	Info<< "Calculating pVelocityDotNormal" << nl << endl;

	surfaceScalarField pVelocityDotNormal
	(
		IOobject
		(
			"pVelocityDotNormal",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
	linearInterpolate(pVelocity) & mesh.Sf()/mag(mesh.Sf())
	);

	Info<< "Calculating ae" << nl << endl;

	surfaceScalarField ae 
	(
		IOobject
		(
			"ae",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		pos(eVelocityDotNormal)
	);

	Info<< "Calculating ap" << nl << endl;

	surfaceScalarField ap 
	(
		IOobject
		(
			"ap",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		pos(pVelocityDotNormal)
	);

	Info<< "Calculating vThermalIons" << nl << endl;

	surfaceScalarField vThermalIons
	(
		IOobject
		(
		"vThermalIons",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
		),
		mesh,
		dimensionedScalar("vThermalIons", dimensionSet(0, 1, -1, 0, 0, 0, 0), Foam::scalar(509.))
	);

	Info << "Calculating vThermalElectrons" << nl << endl;

	surfaceScalarField vThermalElectrons
	(
		IOobject
		(
			"vThermalElectrons",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mesh,
		dimensionedScalar("vThermalElectrons", dimensionSet(0, 1, -1, 0, 0, 0, 0), Foam::scalar(7.08e5))
	);

	Info << "Calculating unitySurfaceField" << nl << endl;

	surfaceScalarField unitySurfaceField
	(
		IOobject
		(
			"unitySurfaceField",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mesh,
		dimensionedScalar("unitySurfaceField", dimensionSet(0, 0, 0, 0, 0, 0, 0), Foam::scalar(1.))
	);

	Info << "Calculating pDensityBoundaryCondition" << nl << endl;

	surfaceScalarField pDensityBoundaryCondition
	(
		IOobject
		(	
			"pDensityBoundaryCondition",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		scalarIonDiffusiveFlux/(0.5*vThermalIons + (2.*ap - 2.*unitySurfaceField)*pVelocityDotNormal)
	);

	Info<< "Calculating gammaIons" << nl << endl;

	volVectorField gammaIons
	(
		IOobject
		(
			"gammaIons",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		pDensity*pVelocity - Dp*fvc::grad(pDensity)
	);

	scalar gammaP = 0.1;
	scalar pi = 3.14159265359;

	dimensionedScalar massIons
	(
		"massIons",
		dimensionSet(1,0,0,0,0,0,0),
		scalar(1.66e-27*29.)
	);

	dimensionedScalar massNeutrals
	(
		"massNeutrals",
		dimensionSet(1,0,0,0,0,0,0),
		scalar(1.66e-27*29.)
	);


	Info<< "Calculating gammaElectronDensity" << nl << endl;

/*	surfaceScalarField gammaElectronDensity
	(
		IOobject
		(
			"gammaElectronDensity",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		(ae - unitySurfaceField) * gammaP * (linearInterpolate(gammaIons) & mesh.Sf()/mag(mesh.Sf())) / eVelocityDotNormal
	); */

	//Dividing by eVelocityDotNormal was a source of a FoamError when the application was run. I believe this could be because we might be dividing by zero

	surfaceScalarField gammaElectronDensity
	(
		IOobject
		(
			"gammaElectronDensity",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		(unitySurfaceField - ae) / eMobility * gammaP * mag( (2.*ap - unitySurfaceField) + 0.5 * sqrt( (8. * (massIons + massNeutrals) * massNeutrals) / (pi * (5.*massIons + 3.*massNeutrals) * massIons) ) ) * pMobility * linearInterpolate(pDensity)
	);

	Info << "Calculating eDensityBoundaryCondition" << nl << endl;

	surfaceScalarField eDensityBoundaryCondition
	(
		IOobject
		(
			"eDensityBoundaryCondition",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		(-scalarElectronDiffusiveFlux - 0.5*vThermalElectrons*gammaElectronDensity - 2.0*(unitySurfaceField - ae)*gammaP*(linearInterpolate(gammaIons) & mesh.Sf()/mag(mesh.Sf())))/((2.0*unitySurfaceField - 2.0*ae) * eVelocityDotNormal - 0.5*vThermalElectrons)
	);
		

//	scalar ae = (average((linearInterpolate(eVelocity) & mesh.Sf())) > test) ? 1. : 0.;
//	scalar ap = (average((linearInterpolate(pVelocity) & mesh.Sf())) > test) ? 1. : 0.;	



/*	surfaceScalarField zeroSurfaceField
	(
		IOobject
		(
			"test",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mesh,
		dimensionedScalar("test", dimensionSet(0, 1, -1, 0, 0, 0, 0), Foam::scalar(0.))
	); */


/*	dimensionedScalar test
	(
		"test",
		dimensionSet(0,3,-1,0,0,0,0),
		scalar(0.0)
	); */	

	Info << "Reached end of createFields.H! " << nl << endl;
		
