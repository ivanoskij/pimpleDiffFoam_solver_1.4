
//////////////////////////////////////////////////////////////////////////////////////
// D-SQUARED LAW /////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

Info<< "Reading field dW\n" << endl;
volScalarField dW	//diameter 
(
    IOobject
    (
        "dW",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh 
);

//////////////////////////////////////////////////////////////////////////////////////

Info<< "Droplets evaporating\n" << endl;

const dimensionedScalar dummySmall1("dummySmall1",dimensionSet(0, 1, -1, 0, 0),1e-5);	//dummy velocity
const dimensionedScalar dummySmall2("dummySmall2",dimensionSet(0, 1, 0, 0, 0),1e-9);	//dummy lenght (1 nm)

//const dimensionedScalar ("dummyT",dimensionSet(0, 0, 1, 0, 0),1);			//dummy time  

//dimensionedScalar tau ("tau", dimensionSet(0, 0, 1, 0, 0), 1);			//time (initialise at zero)
//tau=dummyT*runTime; 	

//////////////////////////////////////////////////////////////////////////////////////

volScalarField Psurf ("Psurf", (Psat*exp(4*sigma*M/((dW+dummySmall2)*R*Tsurf*rhoW)))); 	//vapour pressure at droplet surface (Kelvin equation)

dimensionedScalar Pinf ("Pinf", dimensionSet(1, -1, -2, 0, 0), 0);			//free stream vapour partial pressure (initialise at zero)
Pinf = RH*Psat/100;									//free stream vapour partial pressure

volScalarField Re ("Re", ((max(mag(dummySmall1),mag(U)))*dW/nu));			//Reynolds number
				
const dimensionedScalar Sc ("Sc", nu/Dw);						//Schmidt number

volScalarField Sh ("Sh", (2+0.552*pow(Re,0.5)*pow(Sc,0.333))); 				//Sherwood number

//////////////////////////////////////////////////////////////////////////////////////

volScalarField dW1	//temporary field 1
(
    IOobject
    (
        "dW1",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    2*Dw*M*(Pinf/Tinf-Psurf/Tsurf)*pow(rhoW*R,-1) 
);

volScalarField dW2	//temporary field 2
(
    IOobject
    (
        "dW2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    pow((pow(dW,2)+0.5*dW1*Sh*runTime.time().deltaT().value()),0.5)
);

Info<< "Computing new field dW\n" << endl;
dW += dW2; 

//////////////////////////////////////////////////////////////////////////////////////

/*

dW = pow((pow(dW.oldTime(),2)-Sh*Dw*M*(Psurf/Tsurf-Pinf/Tinf)*pow(rhoW*R,-1)*runTime.value()),0.5);	//d²-law


mathematically: 

	D²(n+1) = D²(n) + something*tau
	hence
	D(n+1) = sqrt [D²(n) + something*tau]

computationally:

	1) compute "something" from previous iteration (n)
	   volScalarField ddW
	   D(n) is known from previous iteration
	   tau is fixed

	2) create new field with the same units of dW
	   volScalarField dddW(n)
	   dddW = sqrt [pow(dW,2) + ddW*tau]

	3) compute new diameter dW
	   dW += dddW

*/
//////////////////////////////////////////////////////////////////////////////////////


/*
int i=0;
int maxIter=2;
double tol=1e-5;
bool ind=false;

dimensionedScalar tau ("tau", dimensionSet(0, 0, 1, 0, 0), 1);
tau=dummyT*runTime; 		

//double tau=runTime.value();

volScalarField magU ("magU", mag(U));

volScalarField ddW("ddW", Sh*Dw*M*(Psurf/Tsurf-Pinf/Tinf)*pow(rhoW*R,-1)*tau);
volScalarField dW_("dW_", dW*0);

dW.storePrevIter();

volScalarField dWprev("dWprev", dW.prevIter());

forAll(mesh.cells(),celli) // or mesh.C() ? 
{
	if(dW[celli]!=0)
	{
		while(ind!=false || i!=maxIter)
		{

		i++;

		Psurf[celli]=Psat[celli]*exp(4.0*sigma.value()*M.value())/(dW[celli]*R.value()*Tsurf.value()*rhoW.value()));
		Re[celli]=max(mag(dummySmall).value(),magU[celli])*dW[celli]/nu.value();
		Sh[celli]=2+0.552*pow(Re[celli],0.5)*pow(Sc,0.333);
		ddW[celli]=Sh*Dw*M*(Psurf[celli]/Tsurf-Pinf/Tinf)*pow(rhoW*R,-1)*tau.value();

			if(dWprev[celli]-ddW[celli])<0)
			{
				ddW[celli]=0
			}

			dW_[celli]==pow(pow(dWprev[celli],2)-ddW,0.5);	//new value
		
			if(mag(dW[celli]-dW_[celli])<=tol)
			{
				ind=true;
				return ind;
			}
			else
			{
				dW[celli]=dW_[celli];	//initialise new value with old value
			}
		
		}
	}
};

Info<< "dW computation completed\n" << endl;
*/

//////////////////////////////////////////////////////////////////////////////////////
/*
// it turns but wrongly
//declarations and initialisation

volScalarField dW_ ("dW_", dW);
volScalarField dW__ ("dW__", dW*0);
volScalarField ddW ("ddW", dW*dW*0);
volScalarField ddW_ ("ddW_", dW*dW*0);
volScalarField delta ("delta", dW*0);

volScalarField Psurf_ ("Psurf_", Psurf);
volScalarField Re_ ("Re_", Re);
volScalarField Sh_ ("Sh_", Sh);

dimensionedScalar tau ("tau", dimensionSet(0, 0, 1, 0, 0), 1);
tau=dummyT*runTime.value();

dW_=dW;
Psurf_=Psurf;
Re_=Re;
Sh_=Sh;

bool ind=false;
int i=0;
int maxIter=2;
double eps=1e-5;

//diameter computation
dW.storePrevIter();
Info<< "computing dW\n" << endl;
while (i!=maxIter || ind)
{
	i++;

	forAll (dW,celli)
	{
		if(dW[celli] == 0)
		{
			Info<< "diameter not found\n" << endl;
			break;
		}
	}
	
	Info<< "up here 1\n" << endl;
	Psurf_=Psat*exp(4*sigma*M/((dW_+dummyBig)*R*Tsurf*rhoW));
	Info<< "up here 2\n" << endl;
	Re_=((max(mag(dummySmall),mag(U)))*dW_/nu);
	Info<< "up here 3\n" << endl;
	Sh_=2+0.552*pow(Re_,0.5)*pow(Sc,0.333);
	Info<< "up here 4\n" << endl;
	ddW=Sh*Dw*M*(Psurf_/Tsurf-Pinf/Tinf)*pow(rhoW*R,-1)*tau;
	Info<< "up here 5\n" << endl;
	ddW_=pow(dW.prevIter(),2) - ddW;
	Info<< "up here 6\n" << endl;
	delta=mag(dW_-dW__);
	Info<< "up here 7\n" << endl;

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
			//Info<< "dW computation converged\n" << endl;
		}
		else
		{
			dW_=dW__;
		}
	}
};
*/
//////////////////////////////////////////////////////////////////////////////////////

//never compiled but it may not work
/*
bool ind=false;

volScalarField ddW ("ddW", Sh*Dw*M*(Psurf/Tsurf-Pinf/Tinf)*pow(rhoW*R,-1)*runTime.value());
volScalarField ddW_ ("ddW_", pow(dW.prevIter(),2) - ddW);

dW.storePrevIter();

fvScalarMatrix dWEqn
{
  fvm::pow(dW,2)
- fvm::pow(dW.prevIter(),2)
==
  fvOptions(dW)
};

forAll(dW, celli)
{
	if (dW[celli] != 0)
	{
		ind=true;
	}
};

forAll(ddW, celli)
{
	if (ddW_[celli] < 0)
	{
		ddW[celli]==0;
	}
};

if(ind)
{
	solve(dWEqn == fvm::ddW);
	fvOptions.correct(dW);
};
*/
//////////////////////////////////////////////////////////////////////////////////////

