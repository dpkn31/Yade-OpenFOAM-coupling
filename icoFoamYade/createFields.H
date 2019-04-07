#include "readGravitationalAcceleration.H"

IOdictionary transportProperties
{
  IOobject
  (
   "transportProperties", 
   runTime.constant(), 
   mesh, 
   IOobject::MUST_READ_IF_MODIFIED, 
   IOobject::NO_WRITE )
}; 


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

dimensionedScalar nu (
 "nu", 
 dimViscosity, 
 transportProperties.lookup("nu") 
); 

dimensionedScalar partDensity (
    "partDensity", 
     dimDensity, 
     transportProperties.lookup("partDensity")
); 

dimensionedScalar fluidDensity( 
    "fluidDensity",
    dimDensity, 
    transportProperties.lookup("fluidDensity")
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


volScalarField alphac
(
    IOobject
    (
      "alpha",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimless, 0)
);


volVectorField divT(
      IOobject(
      "divT",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE),
      fvc::laplacian(U)
);



volTensorField vGrad(
    IOobject(
      "vGrad",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE),
      fvc::grad(U)
    );


volVectorField gradP(
      IOobject(
      "gradP",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE),
      fvc::grad(p)
);



volScalarField uSourceDrag(
    IOobject(
      "uSourceDrag",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("dragCoeff", dimensionSet(0,0,-1,0,0,0,0), 0 )
);


volVectorField uSource(
    IOobject(
      "uSource",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE),
    mesh,
    dimensionedVector("volSource", dimensionSet(0,1,-2,0,0,0,0), vector(0,0,0) )
);



volVectorField uParticle(
    IOobject(
      "uParticle",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE),
    mesh,
    dimensionedVector("uParticle", dimensionSet(0,1,-1,0,0,0,0), vector(0,0,0) )
);


#include "createPhi.H"


volVectorField ddtU_f (
    IOobject (
      "ddtU_f",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE),
      fvc::ddt(U)+fvc::div(phi, U)
);



label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


//singlePhaseTransportModel laminarTransport(U, phi);
//
//autoPtr<incompressible::turbulenceModel> turbulence
//(
//    incompressible::turbulenceModel::New(U, phi, laminarTransport)
//);
//
//#include "createMRF.H"
//#include "createFvOptions.H"
