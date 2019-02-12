
//////////////////////////////////////////////////////////////////////////////////////
// D-SQUARED LAW /////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

volScalarField dW	//diameter 
(
    IOobject
    (
        "dW",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh 
);

//////////////////////////////////////////////////////////////////////////////////////

volScalarField Psurf ("Psurf", (Psat*exp(4*sigma*M/(dW*R*Tsurf*rhoW)))); 		//vapour pressure at droplet surface (Kelvin equation)

dimensionedScalar Pinf ("Pinf", p.dimensions(), 0);		//free stream vapour partial pressure //initialise at zero
Pinf = RH/Psat;							//free stream vapour partial pressure//formula

//////////////////////////////////////////////////////////////////////////////////////
// ADIMENSIONAL NUMBERS //////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

const dimensionedScalar dummyV("dummyV",dimensionSet(0, 1, -1, 0, 0),1e-10);	//minimal velocity allowed for Reloc calculation

volScalarField Reloc ("Reloc", ((max(mag(dummyV),mag(U)))*dW*rhoW/mu));		//local reynolds number
					
const dimensionedScalar Sc ("Sc", nu/Dw);					//Schmidt number

volScalarField Sh ("Sh", (2+0.552*pow(Reloc,0.5)*pow(Sc,0.333))); 		//Sherwood number

//////////////////////////////////////////////////////////////////////////////////////

dW = pow((pow(dW.oldTime(),2)-Sh*Dw*M*(Psurf/Tsurf-Pinf/Tinf)*pow(rhoW*R,-1)*runTime.value()),0.5);	//d²-law

//////////////////////////////////////////////////////////////////////////////////////
