
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

volScalarField Reloc ("Reloc", ((max(mag(dummySmall1),mag(U)))*(dW+dummySmall2)*rhoW/mu));		//local Reynolds number

volScalarField Cw ("Cw", (24/Reloc+4*pow(Reloc,-0.333)));	//drag coefficient//[adim]

volScalarField relt ("relt",(Cw*rhoW*pow(dW,2)/(18*mu)));	//relaxation time//[s]

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Info<< "\n Reading gravitational properties" << endl;		
IOdictionary gravitationalProperties				//gravity
(
    IOobject
    (
        "gravitationalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

const dimensionedVector g(gravitationalProperties.lookup("g"));
const dimensionedScalar magg = mag(g);
const dimensionedVector gHat = g/magg;

volVectorField grav
(
    IOobject
    (
        "grav",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    g
);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Info<< "Computing settling velocity\n" << endl;
volVectorField Uw						//settling velocity
(
    IOobject
    (
        "Uw",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    relt*grav
);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Info<< "computing face flux field phiW\n" << endl;
surfaceScalarField phiW						//settling velocity flux
(
    IOobject
    (
        "phiW",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::flux(Uw)
);
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
