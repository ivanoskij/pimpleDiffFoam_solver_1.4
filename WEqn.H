volScalarField dW3	//temporary field 3
(
    IOobject
    (
        "dW3",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    pow(dW, 3)
);

fvScalarMatrix WEqn
(
    fvm::ddt(W)
  + fvm::div(phi, W)
  + fvm::div(phiW, W)
  - fvm::laplacian(turbulence->nuEff()/Sc,W)
  ==
    fvOptions(W)
);

WEqn.relax();

fvOptions.constrain(WEqn);

solve(WEqn == W*pow(dW, -3)*fvm::ddt(dW3));

fvOptions.correct(W);



