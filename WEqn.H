
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
WEqn.solve();


