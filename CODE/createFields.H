    if (mesh.nCells() != 1)
    {
        FatalErrorInFunction
            << "Solver only applicable to single cell cases"
            << exit(FatalError);
    }

    Info<< "Reading initial conditions.\n" << endl;
    IOdictionary initialConditions
    (
        IOobject
        (
            "initialConditions",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    scalar p0 = readScalar(initialConditions.lookup("p"));
    scalar T0 = readScalar(initialConditions.lookup("Toxidizer"));

    // write base thermo fields - not registered since will be re-read by
    // thermo package

    Info<< "Creating base fields for time " << runTime.timeName() << endl;
    {
        volScalarField Ydefault
        (
            IOobject
            (
                "Ydefault",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar(dimless, 1)
        );

        Ydefault.write();

        volScalarField p
        (
            IOobject
            (
                "p",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar(dimPressure, p0)
        );

        p.write();

        volScalarField T
        (
            IOobject
            (
                "T",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar(dimTemperature, T0)
        );

        T.write();
    }


    Info<< "Reading thermophysical properties\n" << endl;
    autoPtr<rhoReactionThermo> pThermo(rhoReactionThermo::New(mesh));
    rhoReactionThermo& thermo = pThermo();
    thermo.validate(args.executable(), "h");

    basicSpecieMixture& composition = thermo.composition();
    PtrList<volScalarField>& Y = composition.Y();
    volScalarField& p = thermo.p();
    const label nSpecie = Y.size();


    IOField<scalar> Z
    (   
        IOobject
        (
            "Z",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        )
    );

    PtrList<IOField<scalar>> writeData(nSpecie+2);//Y, T, rho
    
    forAll(Y, speciesI)
    {
        writeData.set
        (
            speciesI,
            new IOField<scalar>
            (   
                IOobject
                (
                    Y[speciesI].name(),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                0
            )
        );
    }
    
    writeData.set
    (
        nSpecie,
        new IOField<scalar>
        (   
            IOobject
            (
                "T",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            0
        )
    );
    
    writeData.set
    (
        nSpecie+1,
        new IOField<scalar>
        (   
            IOobject
            (
                "rho",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            0
        )
    );
