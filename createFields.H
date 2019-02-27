
Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//////////////////////////////////////////////////////////////////////////////////////

Info<< "Reading field W\n" << endl;
volScalarField W
(
    IOobject
    (
        "W",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
//////////////////////////////////////////////////////////////////////////////////////

#include "allVariables.H"
#include "evaporation.H"

Info<< "Reading sink term\n" << endl;
volScalarField mW
(
    IOobject
    (
        "mW",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    //mesh,
    rhoW*3.14159/6*pow(dW,3)
);
//////////////////////////////////////////////////////////////////////////////////////

#include "createPhi.H"   

//////////////////////////////////////////////////////////////////////////////////////

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

#include "createMRF.H"
