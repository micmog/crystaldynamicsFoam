#include "createRDeltaT.H"

dimensionedScalar pi("pi", dimensionSet(0, 0, 0, 0, 0), constant::mathematical::pi);

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
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


// volVectorField D
// (
//     IOobject
//     (
//         "D",
//         runTime.timeName(),
//         mesh,
//         IOobject::MUST_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh
// );

const label N_Ori=3;

double phi1[N_Ori]={0.0};

double phi2[N_Ori]={0.0};

double phi3[N_Ori]={0.0};


phi1[0]=0.0;
phi2[0]=0.0;
phi3[0]=0.0;

phi1[1]=0.0;
phi2[1]=0.0;//0.25*pi.value();
phi3[1]=0.25*pi.value();

phi1[2]=0.0;
phi2[2]=0.0;
phi3[2]=0.125*pi.value();




volTensorField omega_dot
(
    IOobject
    (
        "omega_dot",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
        mesh,
dimensionedTensor("omega_dot",dimensionSet(0, 0, -1, 0, 0),tensor::zero)//,
               // zeroGradientFvPatchTensorField::typeName
);

volSymmTensorField epsilon_dot
(
    IOobject
    (
        "epsilon_dot",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
        mesh,
dimensionedSymmTensor("epsilon_dot",dimensionSet(0, 0, -1, 0, 0),symmTensor::zero)//,
               // zeroGradientFvPatchTensorField::typeName
);


// volTensorField gradU
// (
//     IOobject
//     (
//         "gradU",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//         mesh,
// dimensionedTensor("gradU",dimensionSet(0, 0, -1, 0, 0),tensor::zero)//,
//                // zeroGradientFvPatchTensorField::typeName
// );


// volTensorField epsilon
// (
//     IOobject
//     (
//         "epsilon",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//         mesh,
// dimensionedTensor("epsilon",dimensionSet(0, 0, 0, 0, 0),tensor::zero)//,
//                // zeroGradientFvPatchTensorField::typeName
// );

volTensorField s
(
    IOobject
    (
        "s",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
        mesh,
dimensionedTensor("s",dimensionSet(1, -1, -2, 0, 0),tensor::zero)//,
                //zeroGradientFvPatchTensorField::typeName
);

volVectorField divs
(
    IOobject
    (
        "divs",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
        mesh,
dimensionedVector("divs",dimensionSet(1, -2, -2, 0, 0),vector::zero)//,
                //zeroGradientFvPatchTensorField::typeName
);


// volTensorField testSIG
// (
//     IOobject
//     (
//         "testSIG",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//         mesh,
// dimensionedTensor("testSIG",dimensionSet(1, -1, -2, 0, 0),tensor::zero)//,
//                 //zeroGradientFvPatchTensorField::typeName
// );


volTensorField dsdt
(
    IOobject
    (
        "dsdt",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
        mesh,
dimensionedTensor("dsdt",dimensionSet(1, -1, -3, 0, 0),tensor::zero)//,
                //zeroGradientFvPatchTensorField::typeName
);

volTensorField dsdt_crystal
(
    IOobject
    (
        "dsdt_crystal",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
        mesh,
dimensionedTensor("dsdt_crystal",dimensionSet(1, -1, -3, 0, 0),tensor::zero)//,
                //zeroGradientFvPatchTensorField::typeName
);

#include "createPhi.H"


Info<< "Reading phaseProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

autoPtr<twoPhaseChangeModel> phaseChangePtr
(
    twoPhaseChangeModel::New(mixture)
);


// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh
);
// rho.oldTime();

// rho=mixture.rho();

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

pressureReference pressureReference(p, p_rgh, pimple.dict());

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pressureReference.refValue()
      - getRefCellValue(p, pressureReference.refCell())
    );
    p_rgh = p - rho*gh;
}

mesh.schemes().setFluxRequired(p_rgh.name());

volScalarField& alpha1(mixture.alpha1());
mesh.schemes().setFluxRequired(alpha1.name());


#include "createAlphaFluxes.H"

// Construct incompressible turbulence model
incompressibleInterPhaseTransportModel turbulence
(
    U,
    phi,
    alphaPhi1,
    mixture
);


#include "createMRF.H"
#include "createFvModels.H"
#include "createFvConstraints.H"
