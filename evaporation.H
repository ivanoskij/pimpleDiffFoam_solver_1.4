
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

volScalarField Psurf ("Psurf", (Psat*exp(4*sigma*M/(dW*R*Tsurf*rhoW)))); 	//vapour pressure at droplet surface (Kelvin equation)

dimensionedScalar Pinf ("Pinf", p.dimensions(), 0);				//free stream vapour partial pressure //initialise at zero
Pinf = RH/Psat;									//free stream vapour partial pressure//formula

const dimensionedScalar dummyV("dummyV",dimensionSet(0, 1, -1, 0, 0),1e-10);	//minimal velocity allowed for Reloc calculation

volScalarField Re ("Re", ((max(mag(dummyV),mag(U)))*dW/nu));			//Reynolds number
					
const dimensionedScalar Sc ("Sc", nu/Dw);					//Schmidt number

volScalarField Sh ("Sh", (2+0.552*pow(Re,0.5)*pow(Sc,0.333))); 			//Sherwood number

//////////////////////////////////////////////////////////////////////////////////////

//dW = pow((pow(dW.oldTime(),2)-Sh*Dw*M*(Psurf/Tsurf-Pinf/Tinf)*pow(rhoW*R,-1)*runTime.value()),0.5);	//d²-law

//////////////////////////////////////////////////////////////////////////////////////

//declarations and initialisation

volScalarField dW_ ("dW_", dW);
volScalarField dW__ ("dW__", dW*0);
volScalarField ddW ("ddW", dW*0);
volScalarField ddW_ ("ddW_", dW*0);
volScalarField delta ("delta", dW*0);

volScalarField Psurf_ ("Psurf_", Psurf);
volScalarField Re_ ("Re_", Re);
volScalarField Sh_ ("Sh_", Sh);

dW_=dW;
Psurf_=Psurf;
Re_=Re;
Sh_=Sh;

bool ind=false;
int i=0;
int maxIter=1000;
double eps=1e-5;

//diameter computation

while (i!=maxIter || ind!=false)
{
	i++;

	Psurf_=Psat*exp(4*sigma*M/(dW_*R*Tsurf*rhoW));
	Re_=((max(mag(dummyV),mag(U)))*dW_/nu);
	Sh_=2+0.552*pow(Re_,0.5)*pow(Sc,0.333);
	ddW=Sh*Dw*M*(Psurf_/Tsurf-Pinf/Tinf)*pow(rhoW*R,-1)*runTime.value();
	ddW_=pow(dW.prevIter(),2) - ddW;
	delta=mag(dW_-dW__);

	forAll(ddW_, celli)
	{
		if (ddW_[celli] >= 0)
		{
			dW__=pow(pow(dW.prevIter(),2)-ddW,0.5);
		}
		else
		{
			dW__=pow(pow(dW.prevIter(),2),0.5);
		}
	}
	
	forAll(delta, celli)
	{
		if (delta[celli] <= eps)
		{
			ind=true;
			dW=dW__;
			Info<< "converged\n" << endl;
		}
		else
		{
			dW_=dW__;
		}
	}
};

//////////////////////////////////////////////////////////////////////////////////////

