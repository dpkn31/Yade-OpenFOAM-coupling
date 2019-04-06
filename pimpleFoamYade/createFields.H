#include "readGravitationalAcceleration.H"

word continuousPhaseName
(
    IOdictionary
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ
        )
    ).lookup("continuousPhaseName")
);


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




Info<< "Reading field U\n" << endl;
volVectorField Uc
(
    IOobject
    (
        IOobject::groupName("U", continuousPhaseName),
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

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

dimensionedScalar partDensity (
    "partDensity",
     dimDensity,
     transportProperties.lookup("partDensity")
);


Info<< "Reading/calculating continuous-phase face flux field phic\n"
    << endl;

surfaceScalarField phic
(
    IOobject
    (
        IOobject::groupName("phi", continuousPhaseName),
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(Uc) & mesh.Sf()
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());

Info<< "Creating turbulence model\n" << endl;

singlePhaseTransportModel continuousPhaseTransport(Uc, phic);

dimensionedScalar rhocValue
(
    IOobject::groupName("rho", continuousPhaseName),
    dimDensity,
    continuousPhaseTransport.lookup
    (
        IOobject::groupName("rho", continuousPhaseName)
    )
);


dimensionedScalar nuValue (
   "kinematic viscosity",
  dimViscosity,
  transportProperties.lookup("nu")
);

volScalarField rhoc
(
    IOobject
    (
        rhocValue.name(),
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    rhocValue
);

volScalarField muc
(
    IOobject
    (
        IOobject::groupName("mu", continuousPhaseName),
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    rhoc*continuousPhaseTransport.nu()
);

Info << "Creating field alphac\n" << endl;
// alphac must be constructed before the cloud
// so that the drag-models can find it
volScalarField alphac
(
    IOobject
    (
        IOobject::groupName("alpha", continuousPhaseName),
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
      fvc::laplacian(Uc)
);



volTensorField vGrad(
    IOobject(
      "vGrad",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE),
      fvc::grad(Uc)
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


volVectorField ddtU_f (
    IOobject (
      "ddtU_f",
      runTime.timeName(),
      mesh,
      IOobject::NO_READ,
      IOobject::NO_WRITE),
      fvc::ddt(Uc)+fvc::div(phic, Uc)
);



surfaceScalarField alphacf("alphacf", fvc::interpolate(alphac));

surfaceScalarField alphaPhic
(
    IOobject::groupName
    (
        "alphaPhi",
        continuousPhaseName
    ),
    alphacf*phic
);

autoPtr<PhaseIncompressibleTurbulenceModel<singlePhaseTransportModel>>
continuousPhaseTurbulence
(
    PhaseIncompressibleTurbulenceModel<singlePhaseTransportModel>::New
    (
        alphac,
        Uc,
        alphaPhic,
        phic,
        continuousPhaseTransport
    )
);
