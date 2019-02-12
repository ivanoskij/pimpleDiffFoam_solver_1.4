///////////////////////////////////////////////////////////////////////////
// DROPLETS AND WATER /////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
dimensionedScalar Dw		//diffusivity of water vapour in air
(
    "Dw",
    transportProperties.lookup("Dw")
);
dimensionedScalar M		//molecular weight of water
(
    "M",
    transportProperties.lookup("M")
);
dimensionedScalar R		//gas coostant 
(
    "R",
    transportProperties.lookup("R")
);

dimensionedScalar Tsurf		//temperature at droplet surface
(
    "Tsurf",
    transportProperties.lookup("Tsurf")
);

dimensionedScalar sigma		//surface tension
(
    "sigma",		
    transportProperties.lookup("sigma")
);
dimensionedScalar rhoW 		//density of water
(
    "rhoW",
    transportProperties.lookup("rhoW")
);
dimensionedScalar Psat		//vapour saturation pressure 
(
    "Psat",
    transportProperties.lookup("Psat")
);
///////////////////////////////////////////////////////////////////////////
// AIR ////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

dimensionedScalar mu		//air viscosity
(
    "mu",
    transportProperties.lookup("mu")
);

dimensionedScalar nu		//air kinematic viscosity
(
    "nu",
    transportProperties.lookup("nu")
);

dimensionedScalar RH		//relative humidity
(
    "RH",
    transportProperties.lookup("RH")
);

dimensionedScalar Tinf		//free stream temperature 
(
    "Tinf",
    transportProperties.lookup("Tinf")
);
///////////////////////////////////////////////////////////////////////////
