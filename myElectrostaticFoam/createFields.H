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

    dimensionedScalar epsilon0
    (
        physicalProperties.lookup("epsilon0")
    );

    dimensionedScalar muip
    (
        physicalProperties.lookup("muip")
    );

    dimensionedScalar muem
    (
        physicalProperties.lookup("muem")
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


    Info<< "Reading field ip\n" << endl;
    volScalarField ip
    (
        IOobject
        (
            "ip",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field em\n" << endl;
    volScalarField em
    (
        IOobject
        (
            "em",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Calculating field ipFlux\n" << endl;
    surfaceScalarField ipFlux
    (
        IOobject
        (
            "ipFlux",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        -muip*mesh.magSf()*fvc::snGrad(phi)
    );

    Info<< "Calculating field emFlux\n" << endl;
    surfaceScalarField emFlux
    (
        IOobject
        (
            "emFlux",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        muem*mesh.magSf()*fvc::snGrad(phi)
    );
