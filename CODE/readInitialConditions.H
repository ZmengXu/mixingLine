/*
    scalar minZ=0.0;
    
    scalar maxZ=1.0;

    initialConditions.lookup("min:maxZ") >> minZ >> maxZ;

    forAll(Z, celli)
    {
        Z[celli] = minZ + (maxZ - minZ)*celli/Z.size();
    }
*/
    scalar Tfuel = readScalar(initialConditions.lookup("Tfuel"));

    scalar Toxidizer = readScalar(initialConditions.lookup("Toxidizer"));
    
    scalarList Yfuel(nSpecie, 0.0);
    
    dictionary fuelDict(initialConditions.subDict("fuel"));
    
    forAll(Y, i)
    {
        const word& name = Y[i].name();
        if (fuelDict.found(name))
        {
            Yfuel[i] = readScalar(fuelDict.lookup(name));
        }
    }
    const scalar mTotFuel = sum(Yfuel);
    forAll(Y, i)
    {
        Yfuel[i] /= mTotFuel;
    }
    
    
    scalarList Yoxidizer(nSpecie, 0.0);
    
    dictionary oxidizerDict(initialConditions.subDict("oxidizer"));
    
    forAll(Y, i)
    {
        const word& name = Y[i].name();
        if (oxidizerDict.found(name))
        {
            Yoxidizer[i] = readScalar(oxidizerDict.lookup(name));
        }
    }
    const scalar mTotOxidizer = sum(Yoxidizer);
    forAll(Y, i)
    {
        Yoxidizer[i] /= mTotOxidizer;
    }
    
    const scalar hfuel = ::h0(thermo, Yfuel, p[0], Tfuel);

    const scalar hoxidizer = ::h0(thermo, Yoxidizer, p[0], Toxidizer);
    
