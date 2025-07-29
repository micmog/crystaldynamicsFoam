Info << "Reading physicalProperties\n"
     << endl;

IOdictionary physicalProperties(
    IOobject(
        "physicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE));

IOdictionary solutionDict(
    IOobject(
        "fvSolution",
        runTime.system(), // folder where the fvSolution file is located
        runTime,
        IOobject::MUST_READ));

IOdictionary controlDict(
    IOobject(
        "controlDict",
        runTime.system(), // folder where the fvSolution file is located
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE));

dictionary &PimpleDICT = solutionDict.subDict("PIMPLE");

int n_sub_steps(readLabel(PimpleDICT.lookup("nOuterCorrectors")));

// scalar maxTimeStep(readScalar(controlDict.lookup("maxDeltaT")));

// Info<<"TREST: "<<n_sub_steps<<endl;

// exit(0);

dimensionedScalar pi("pi", dimensionSet(0, 0, 0, 0, 0), constant::mathematical::pi);

dimensionedScalar Cv("Cv", dimensionSet(0, 2, -2, -1, 0), physicalProperties.lookup("Cv"));
dimensionedScalar kappa("kappa", dimensionSet(1, 1, -3, -1, 0), physicalProperties.lookup("kappa"));
dimensionedScalar rho_0("rho_0", dimensionSet(1, -3, 0, 0, 0), physicalProperties.lookup("rho_0"));

//  Info<<"TEST"<<(n_alpha[0] * m_alpha[0])<<endl;

// exit(0);

volTensorField D_P(
    IOobject(
        "D_P",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedTensor("D_P", dimensionSet(0, 0, -1, 0, 0), tensor::zero)//,
    // zeroGradientFvPatchTensorField::typeName
    );

volTensorField Omega_P
(
    IOobject
    (
        "Omega_P",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
        mesh,
dimensionedTensor("Omega_P",dimensionSet(0, 0, -1, 0, 0),tensor::zero),
                zeroGradientFvPatchTensorField::typeName
);

// volTensorField Omega_E
// (
//     IOobject
//     (
//         "Omega_E",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//         mesh,
// dimensionedTensor("Omega_E",dimensionSet(0, 0, -1, 0, 0),tensor::zero),
//                 zeroGradientFvPatchTensorField::typeName
// );

// for ()
// {

// }

// const Switch Hardening(physicalProperties.lookupOrDefault<Switch>("Hardening", false));

// // if(Hardening){

// scalar K_0(physicalProperties.lookupOrDefault<scalar>("K_0", 0.0));
// scalar K_n(physicalProperties.lookupOrDefault<scalar>("K_n", 0.0));
// scalar K_B(physicalProperties.lookupOrDefault<scalar>("K_B", 0.0));
// scalar K_Beta_0(physicalProperties.lookupOrDefault<scalar>("K_Beta_0", 0.0));
// scalar K_Beta_1(physicalProperties.lookupOrDefault<scalar>("K_Beta_1", 0.0));

// // dimensionedScalar K_0("K_0",dimensionSet(1, -1, -2, 0, 0),physicalProperties.lookup("K_0"));
// // dimensionedScalar K_n("K_n",dimensionSet(0, 0, 0, 0, 0),physicalProperties.lookup("K_n"));
// // dimensionedScalar K_B("K_B",dimensionSet(0, 0, 0, 0, 0),physicalProperties.lookup("K_B"));
// // dimensionedScalar K_Beta_0("K_Beta_0",dimensionSet(0, 0, 0, 0, 0),physicalProperties.lookup("K_Beta_0"));
// // dimensionedScalar K_Beta_1("K_Beta_1",dimensionSet(0, 0, 0, 0, 0),physicalProperties.lookup("K_Beta_1"));

// // }

// dimensionedScalar sigma_Y_0("sigma_Y_0",dimensionSet(1, -1, -2, 0, 0),physicalProperties.lookup("sigma_Y_0"));

Info << "Reading field U\n"
     << endl;
volVectorField U(
    IOobject(
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE),
    mesh);

volVectorField D
(
    IOobject
    (
        "D",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("D",dimensionSet(0, 1, 0, 0, 0),vector::zero)
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
// dimensionedTensor("gradU",dimensionSet(0, 0, -1, 0, 0),tensor::zero),
//                zeroGradientFvPatchTensorField::typeName
// );

// volTensorField L
// (
//     IOobject
//     (
//         "L",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//         mesh,
// dimensionedTensor("L",dimensionSet(0, 0, -1, 0, 0),tensor::zero),
//                zeroGradientFvPatchTensorField::typeName
// );

// volTensorField gradUT
// (
//     IOobject
//     (
//         "gradUT",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//         mesh,
// dimensionedTensor("gradU",dimensionSet(0, 0, -1, 0, 0),tensor::zero),
//                zeroGradientFvPatchTensorField::typeName
// );




Info << "Reading field Temperature\n"
     << endl;
volScalarField T(
    IOobject(
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE),
    mesh);

volScalarField K(
    IOobject(
        "K",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("K", dimensionSet(0, 2, -2, 0, 0), 0.0));

// volScalarField trsig
// (
//     IOobject
//     (
//         "trsig",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("trsig",dimensionSet(1, -1, -3, 0, 0),0.0)
// );

// volScalarField tredot
// (
//     IOobject
//     (
//         "tredot",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("tredot",dimensionSet(0, 0, -1, 0, 0),0.0)
// );

// volVectorField TESTADV
// (
//     IOobject
//     (
//         "TESTADV",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedVector("TESTADV",dimensionSet(1, -2, -2, 0, 0),vector::zero),
//                 zeroGradientFvPatchVectorField::typeName
// );

// volScalarField DP
// (
//     IOobject
//     (
//         "DP",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("DP",dimensionSet(1, -1, -3, 0, 0),0.0),
//                 zeroGradientFvPatchScalarField::typeName
// );

// volScalarField DP2
// (
//     IOobject
//     (
//         "DP2",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("DP2",dimensionSet(1, -1, -3, 0, 0),0.0),
//                 zeroGradientFvPatchScalarField::typeName
// );

volScalarField rho(
    IOobject(
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("rho", dimensionSet(1, -3, 0, 0, 0), rho_0.value()),
    zeroGradientFvPatchScalarField::typeName);

volScalarField p(
    IOobject(
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("p", dimensionSet(1, -1, -2, 0, 0), 0.0),
    zeroGradientFvPatchScalarField::typeName);

surfaceScalarField rhoPhi(
    IOobject(
        "rhoPhi",
        runTime.timeName(),
        mesh),
    mesh,
    dimensionedScalar(dimMass / dimTime, 0));

// volScalarField divU
// (
//     IOobject
//     (
//         "divU",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::AUTO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("divU",dimensionSet(0, 0, -1, 0, 0),0.0),
//                 zeroGradientFvPatchScalarField::typeName
// );

volScalarField Sum_Phases(
    IOobject(
        "Sum_Phases",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("Sum_Phases", dimensionSet(0, 0, 0, 0, 0), 0.0));

// const label N_Ori=20;

IOdictionary GrainProperties(
    IOobject(
        "GrainProperties",   // dictionary name
        runTime.constant(),  // dict is found in "constant"
        mesh,                // registry for the dict
        IOobject::MUST_READ, // must exist, otherwise failure
        IOobject::NO_WRITE   // dict is only read by the solver
        ));

label N_Ori(readLabel(GrainProperties.lookup("N_Phases")));
Info << "Number of Phase Field Order Parameters: " << N_Ori << endl;

// Ready to read in multiple stiffness matrices

label N_Mats(readLabel(GrainProperties.lookup("N_Materials")));
Info << "Number of Materials: " << N_Mats << endl;

labelList Phase2Material(N_Ori, 0);

GrainProperties.lookup("Phase2Material") >> Phase2Material;

for (label c = 0; c < N_Ori; c++)
{
    Info << "Phase[" << c << "] = Material[" << Phase2Material[c] << "]" << endl;
}

Fastor::Tensor<double, 3, 3, 3, 3> Cprime2[N_Mats]; // array of tensor objects



std::vector<Material> MatList(N_Mats, Material());

for (label i = 0; i < N_Mats; i++)
{

    word Mat_dict("Material." + name(i));
    dictionary Mat = GrainProperties.subDict(Mat_dict);

    Fastor::Tensor<double, 3, 3, 3, 3> Cprime_temp_Read(0.0);

    for (int m = 0; m < 3; m++)
    {
        for (int n = 0; n < 3; n++)
        {
            for (int p = 0; p < 3; p++)
            {
                for (int q = 0; q < 3; q++)
                {

                    word Cijkl("C" + name(m) + name(n) + name(p) + name(q));

                    Cprime_temp_Read(m, n, p, q) = Mat.lookupOrDefault<scalar>(Cijkl, 0.0);
                }
            }
        }
    }

    Cprime2[i] = Cprime_temp_Read;

    label N_slips_Temp(readLabel(Mat.lookup("N_slips")));

    MatList[i].N_Slips = N_slips_Temp;

    for (int j = 0; j < MatList[i].N_Slips; j++)
    {

        word Slip_direction_name = "m" + name(j);
        word Slip_normal_name = "n" + name(j);

        vector temp_m(Mat.lookup(Slip_direction_name));
        vector temp_n(Mat.lookup(Slip_normal_name));

        //   m_alpha[i][j] = temp_m;
        //   n_alpha[i][j] = temp_n;

        MatList[i].SP_norm.push_back(temp_n);
        MatList[i].SP_direction.push_back(temp_m);
    }

    scalar initialCRSS_temp(readScalar(Mat.lookup("initialCRSS")));
    MatList[i].initial_CRSS = initialCRSS_temp;

    scalar tau_S_temp(readScalar(Mat.lookup("tau_S")));
    MatList[i].tau_S = tau_S_temp;

    scalar h_0_temp(readScalar(Mat.lookup("h_0")));
    MatList[i].h_0 = h_0_temp;

    scalar r_temp(readScalar(Mat.lookup("r")));
    MatList[i].r = r_temp;

    scalar k_temp(readScalar(Mat.lookup("k_exponent")));
    MatList[i].k_exponent = k_temp;

    scalar gammadot_0_temp(readScalar(Mat.lookup("gammadot_0")));
    MatList[i].ref_shear_rate = gammadot_0_temp;
}

// Info<<"TEST: "<<MatList[0].SP_norm[0]<<endl;

// Info<<"TEST: "<<MatList[0].SP_norm[11]<<endl;

// exit(0);

PtrList<PtrList<volScalarField>> g_alphas(N_Ori);
forAll(g_alphas, i)
{
    g_alphas.set(i, new PtrList<volScalarField>(MatList[Phase2Material[i]].N_Slips));
    forAll(g_alphas[i], j)
    {
        g_alphas[i].set(
            j,
            new volScalarField(
                IOobject(
                    "g_alpha." + Foam::name(i) + "." + Foam::name(j),
                    runTime.timeName(),
                    mesh,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE),
                mesh,
                dimensionedScalar("g_alphas", dimensionSet(0, 0, 0, 0, 0, 0, 0), MatList[Phase2Material[i]].initial_CRSS /*16.0e6*/), // should be initial for particular material
                zeroGradientFvPatchScalarField::typeName));
    }
}

// PtrList<volScalarField> g_alphas(N_Slips);

// for(label count = 0; count < N_Slips; count++)	//hardening field per slip system
// {

//     word name_g_alpha ("g_alpha." + name(count));

//     g_alphas.set
//     (
//         count,
//         new volScalarField
//         (
//             IOobject
//             (
//                 name_g_alpha,
//                 runTime.timeName(),
//                 mesh,
//                 IOobject::READ_IF_PRESENT,//IOobject::MUST_READ,
//                 IOobject::AUTO_WRITE
//             ),
//             mesh,
//             dimensionedScalar("g_alpha",dimensionSet(0, 0, 0, 0, 0),16.0e6),
//              zeroGradientFvPatchScalarField::typeName
//         )
//     );
// }

// PtrList<volScalarField> g_alpha_dots(N_Slips);//source terms on RHS of G_alpha_equation

// for(label count = 0; count < N_Slips; count++)	//hardening field per slip system
// {

//     word name_g_alpha_dot ("g_alpha_dot." + name(count));

//     g_alpha_dots.set
//     (
//         count,
//         new volScalarField
//         (
//             IOobject
//             (
//                 name_g_alpha_dot,
//                 runTime.timeName(),
//                 mesh,
//                 IOobject::READ_IF_PRESENT,//IOobject::MUST_READ,
//                 IOobject::NO_WRITE
//             ),
//             mesh,
//             dimensionedScalar("g_alpha_dot",dimensionSet(0, 0, -1, 0, 0),0.0),
//              zeroGradientFvPatchScalarField::typeName
//         )
//     );
// }

PtrList<PtrList<volScalarField>> g_alpha_dots(N_Ori);
forAll(g_alpha_dots, i)
{
    g_alpha_dots.set(i, new PtrList<volScalarField>(MatList[Phase2Material[i]].N_Slips));
    forAll(g_alpha_dots[i], j)
    {
        g_alpha_dots[i].set(
            j,
            new volScalarField(
                IOobject(
                    "g_alpha_dots." + Foam::name(i) + "." + Foam::name(j),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE),
                mesh,
                dimensionedScalar("g_alpha_dots", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.0),
                zeroGradientFvPatchScalarField::typeName));
    }
}

// PtrList<volScalarField> gamma_alpha_dots(N_Slips);//source terms on RHS of G_alpha_equation

// for(label count = 0; count < N_Slips; count++)	//hardening field per slip system
// {

//     word name_gamma_alpha_dot ("gamma_alpha_dot." + name(count));

//     gamma_alpha_dots.set
//     (
//         count,
//         new volScalarField
//         (
//             IOobject
//             (
//                 name_gamma_alpha_dot,
//                 runTime.timeName(),
//                 mesh,
//                 IOobject::READ_IF_PRESENT,//IOobject::MUST_READ,
//                 IOobject::AUTO_WRITE
//             ),
//             mesh,
//             dimensionedScalar("gamma_alpha_dot",dimensionSet(0, 0, -1, 0, 0),0.0),
//              zeroGradientFvPatchScalarField::typeName
//         )
//     );
// }

PtrList<PtrList<volScalarField>> gamma_alpha_dots(N_Ori);
forAll(gamma_alpha_dots, i)
{
    gamma_alpha_dots.set(i, new PtrList<volScalarField>(MatList[Phase2Material[i]].N_Slips));
    forAll(gamma_alpha_dots[i], j)
    {
        gamma_alpha_dots[i].set(
            j,
            new volScalarField(
                IOobject(
                    "gamma_alpha_dots." + Foam::name(i) + "." + Foam::name(j),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE),
                mesh,
                dimensionedScalar("gamma_alpha_dots", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.0),
                zeroGradientFvPatchScalarField::typeName));
    }
}

// PtrList<PtrList<volScalarField>> gamma_alphas(N_Ori);
// forAll(gamma_alphas, i)
// {
//     gamma_alphas.set(i, new PtrList<volScalarField>(MatList[Phase2Material[i]].N_Slips));
//     forAll(gamma_alphas[i], j)
//     {
//         gamma_alphas[i].set(
//             j,
//             new volScalarField(
//                 IOobject(
//                     "gamma_alphas." + Foam::name(i) + "." + Foam::name(j),
//                     runTime.timeName(),
//                     mesh,
//                     IOobject::READ_IF_PRESENT,
//                     IOobject::AUTO_WRITE),
//                 mesh,
//                 dimensionedScalar("gamma_alphas", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0),
//                 zeroGradientFvPatchScalarField::typeName));
//     }
// }

volScalarField sumgamma
(
    IOobject
    (
        "sumgamma",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("sumgamma",dimensionSet(0, 0, -1, 0, 0),0.0),
                zeroGradientFvPatchScalarField::typeName
);


// PtrList<volScalarField> total_slip_phase(N_Ori); // pointer to phase fields

// for (label count = 0; count < N_Ori; count++) // how to assign the popbal array of fields
// {

//     word name_ni("total_slip." + name(count));

//     total_slip_phase.set(
//         count,
//         new volScalarField(
//             IOobject(
//                 name_ni,
//                 runTime.timeName(),
//                 mesh,
//                 IOobject::READ_IF_PRESENT, // IOobject::MUST_READ,
//                 IOobject::AUTO_WRITE),
//             mesh,
//             dimensionedScalar("ni", dimensionSet(0, 0, 0, 0, 0), 0.0),
//             zeroGradientFvPatchScalarField::typeName));
// }




PtrList<volScalarField> PopBal(N_Ori); // pointer to phase fields

for (label count = 0; count < N_Ori; count++) // how to assign the popbal array of fields
{

    word name_ni("n." + name(count));

    PopBal.set(
        count,
        new volScalarField(
            IOobject(
                name_ni,
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT, // IOobject::MUST_READ,
                IOobject::AUTO_WRITE),
            mesh,
            dimensionedScalar("ni", dimensionSet(0, 0, 0, 0, 0), 0.0),
            zeroGradientFvPatchScalarField::typeName));
}

PtrList<volScalarField> PopBal_av(N_Ori);

for(label count = 0; count < N_Ori; count++)
{

    word name_ni_av ("n_av" + name(count));

    PopBal_av.set
    (
        count,
        new volScalarField
        (
            IOobject
            (
                name_ni_av,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("niav",dimensionSet(0, 0, 0, 0, 0),0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}

Sum_Phases *= 0.0;
for (label count = 0; count < N_Ori; count++)
{

    volScalarField &ni = PopBal[count];

    Sum_Phases += ni;
}
Sum_Phases.correctBoundaryConditions();

//  Project any out of bounds phases back onto the Gibbs Simplex
for (label count = 0; count < N_Ori; count++)
{

    volScalarField &ni = PopBal[count];

    ni = max(min(ni, 1.0), 0.0);
    volScalarField temp_holder = ni;

    forAll(mesh.C(), celli)
    {
        if (Sum_Phases[celli] > 0.0)
        {
            temp_holder[celli] *= (1.0 / Sum_Phases[celli]);
        }
    }

    forAll(Sum_Phases.boundaryField(), patchI)
{
    if
    (
        !Sum_Phases.boundaryField()[patchI].coupled()
     && Sum_Phases.boundaryField()[patchI].size() > 0
    )
    {
        forAll(Sum_Phases.boundaryField()[patchI], faceI)
        {

            if (Sum_Phases.boundaryField()[patchI][faceI] > 0.0)
        {
            temp_holder.boundaryFieldRef()[patchI][faceI] *= (1.0 / Sum_Phases.boundaryField()[patchI][faceI]);
        }

        }

    }

}



    ni = temp_holder;
}



volScalarField GBREF
(
    IOobject
    (
        "GBREF",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("GBREF",dimensionSet(0,0,0,0,0),0.0)
);



// PtrList<volTensorField> Omega_P(N_Ori); // pointer to phase fields

// for (label count = 0; count < N_Ori; count++) // how to assign the popbal array of fields
// {

//     word name_Omega_P("nOmega_P." + name(count));

//     Omega_P.set(
//         count,
//         new volTensorField(
//             IOobject(
//                 name_Omega_P,
//                 runTime.timeName(),
//                 mesh,
//                 IOobject::READ_IF_PRESENT, // IOobject::MUST_READ,
//                 IOobject::AUTO_WRITE),
//             mesh,
//             dimensionedTensor("Omega_P", dimensionSet(0, 0, -1, 0, 0), tensor::zero)//,
//             // zeroGradientFvPatchTensorField::typeName
//             )
//             );
// }

// Random randObj2(clock::getTime());

// double phi1[N_Ori] = {0.0};

// double phi2[N_Ori] = {0.0};

// double phi3[N_Ori] = {0.0};

List<double> phi1 = GrainProperties.lookup("phi1");
List<double> phi2 = GrainProperties.lookup("phi2");
List<double> phi3 = GrainProperties.lookup("phi3");

// for (label count_j = 0; count_j < N_Ori; count_j++)
// {

//     phi1[count_j] =  randObj2.scalar01()*2.0*pi.value();
//     phi2[count_j] = 0.0; // randObj2.scalar01()*2.0*pi.value();
//     phi3[count_j] = 0.0; // randObj2.scalar01()*2.0*pi.value();
// }

// phi1[0]=0.0;
// phi2[0]=0.0;
// phi3[0]=0.0;

// phi1[1]=2.0*pi.value();
// phi2[1]=0.0;
// phi3[1]=0.0;

// phi1[2]=2.0*pi.value();
// phi2[2]=0.0;
// phi3[2]=0.0;

// phi1[0] = 0.25 * pi.value();
// phi2[0] = 0.0;
// phi3[0] = 0.0;

// PtrList<volScalarField> Rot1(N_Ori); // pointer to phase fields

// for (label count = 0; count < N_Ori; count++) // how to assign the popbal array of fields
// {

//     word name_ROT("ROT.1." + name(count));

//     Rot1.set(
//         count,
//         new volScalarField(
//             IOobject(
//                 name_ROT,
//                 runTime.timeName(),
//                 mesh,
//                 IOobject::READ_IF_PRESENT, // IOobject::MUST_READ,
//                 IOobject::AUTO_WRITE),
//             mesh,
//             dimensionedScalar("ROT.1", dimensionSet(0, 0, 0, 0, 0), phi1[count]),
//             zeroGradientFvPatchScalarField::typeName));
// }

// PtrList<volScalarField> Rot2(N_Ori); // pointer to phase fields

// for (label count = 0; count < N_Ori; count++) // how to assign the popbal array of fields
// {

//     word name_ROT("ROT.2." + name(count));

//     Rot2.set(
//         count,
//         new volScalarField(
//             IOobject(
//                 name_ROT,
//                 runTime.timeName(),
//                 mesh,
//                 IOobject::READ_IF_PRESENT, // IOobject::MUST_READ,
//                 IOobject::AUTO_WRITE),
//             mesh,
//             dimensionedScalar("ROT.2", dimensionSet(0, 0, 0, 0, 0), phi2[count]),
//             zeroGradientFvPatchScalarField::typeName));
// }

// PtrList<volScalarField> Rot3(N_Ori); // pointer to phase fields

// for (label count = 0; count < N_Ori; count++) // how to assign the popbal array of fields
// {

//     word name_ROT("ROT.3." + name(count));

//     Rot3.set(
//         count,
//         new volScalarField(
//             IOobject(
//                 name_ROT,
//                 runTime.timeName(),
//                 mesh,
//                 IOobject::READ_IF_PRESENT, // IOobject::MUST_READ,
//                 IOobject::AUTO_WRITE),
//             mesh,
//             dimensionedScalar("ROT.3", dimensionSet(0, 0, 0, 0, 0), phi3[count]),
//             zeroGradientFvPatchScalarField::typeName));
// }


volTensorField RotationMatrixTemp(
    IOobject(
        "RotationMatrixTemp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    mesh,
    dimensionedTensor("RotationMatrixTemp", dimensionSet(0, 0, 0, 0, 0), tensor::zero),
     zeroGradientFvPatchTensorField::typeName
    );





for (label i = 0; i < N_Ori; i++)
{


 tensor temprot(
    Foam::cos(phi2[i])*Foam::cos(phi3[i]), 
    (Foam::sin(phi1[i])*Foam::sin(phi2[i])*Foam::cos(phi3[i]))-(Foam::cos(phi1[i])*Foam::sin(phi3[i])),//-Foam::cos(phi3)*Foam::sin(phi2), 
    (Foam::cos(phi1[i])*Foam::sin(phi2[i])*Foam::cos(phi3[i]))+(Foam::sin(phi1[i])*Foam::sin(phi3[i])),//Foam::sin(phi2)*Foam::sin(phi3),
   
    Foam::cos(phi2[i])*Foam::sin(phi3[i]), 
    (Foam::sin(phi1[i])*Foam::sin(phi2[i])*Foam::sin(phi3[i]))+(Foam::cos(phi1[i])*Foam::cos(phi3[i])), 
    (Foam::cos(phi1[i])*Foam::sin(phi2[i])*Foam::sin(phi3[i]))-(Foam::sin(phi1[i])*Foam::cos(phi3[i])),
    
    -Foam::sin(phi2[i]), 
    Foam::sin(phi1[i])*Foam::cos(phi2[i]), 
    Foam::cos(phi1[i])*Foam::cos(phi2[i])
);//https://en.wikipedia.org/wiki/Rotation_matrix - 
//"represents an extrinsic rotation whose (improper) Euler angles are α, β, γ, about axes x, y, z" - phi1 ==alpha, phi2==beta,phi3==gamma

RotationMatrixTemp+=PopBal[i]*temprot;

}


volTensorField RotationMatrix(
    IOobject(
        "RotationMatrix",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
        RotationMatrixTemp
    // mesh,
    // dimensionedTensor("RotationMatrix", dimensionSet(0, 0, 0, 0, 0), RotationMatrixTemp)//,
    // zeroGradientFvPatchTensorField::typeName
    );



volScalarField Rot1(
    IOobject(
        "Rot1",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("Rot1", dimensionSet(0, 0, 0, 0, 0), 0.0));

volScalarField Rot2(
    IOobject(
        "Rot2",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("Rot2", dimensionSet(0, 0, 0, 0, 0), 0.0));

    volScalarField Rot3(
    IOobject(
        "Rot3",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("Rot3", dimensionSet(0, 0, 0, 0, 0), 0.0));

// Material TEST();

// Info<<MatList[0].RK4()<<endl;

// exit(0);

volTensorField testtraceeps
(
    IOobject
    (
        "testtraceeps",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
        mesh,
dimensionedTensor("testtraceeps",dimensionSet(1, -1, -3, 0, 0),tensor::zero)//,
            //    zeroGradientFvPatchTensorField::typeName
);

volTensorField omega_dot(
    IOobject(
        "omega_dot",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedTensor("omega_dot", dimensionSet(0, 0, -1, 0, 0), tensor::zero)//,
    // zeroGradientFvPatchTensorField::typeName
    );

volTensorField epsilon_dot(
    IOobject(
        "epsilon_dot",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedTensor("epsilon_dot", dimensionSet(0, 0, -1, 0, 0), tensor::zero)//,
    // zeroGradientFvPatchTensorField::typeName
    );

// volTensorField omega_tot
// (
//     IOobject
//     (
//         "omega_tot",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//         mesh,
// dimensionedTensor("omega_tot",dimensionSet(0, 0, 0, 0, 0),tensor::zero),
//                zeroGradientFvPatchTensorField::typeName
// );

volScalarField Acc_Strain(
    IOobject(
        "Acc_Strain",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("Acc_Strain", dimensionSet(0, 0, 0, 0, 0), 0.0),
    zeroGradientFvPatchScalarField::typeName
    );

volTensorField sigma(
    IOobject(
        "sigma",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedTensor("sigma", dimensionSet(1, -1, -2, 0, 0), tensor::zero)//,
    // zeroGradientFvPatchTensorField::typeName
    );

// volTensorField sigmadiff(
//     IOobject(
//         "sigmadiff",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE),
//     mesh,
//     dimensionedTensor("sisigmadiffgmadiff", dimensionSet(1, -1, -2, 0, 0), tensor::zero)//,
//     // zeroGradientFvPatchTensorField::typeName
//     );

    // volTensorField sigmatemp(
    // IOobject(
    //     "sigmatemp",
    //     runTime.timeName(),
    //     mesh,
    //     IOobject::READ_IF_PRESENT,
    //     IOobject::AUTO_WRITE),
    // mesh,
    // dimensionedTensor("sigmatemp", dimensionSet(1, -1, -2, 0, 0), tensor::zero)//,
    // // zeroGradientFvPatchTensorField::typeName
    // );



volVectorField divsigma(
    IOobject(
        "divsigma",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE),
    mesh,
    dimensionedVector("divsigma", dimensionSet(1, -2, -2, 0, 0), vector::zero) //,
                                                                               // zeroGradientFvPatchTensorField::typeName
);

volVectorField Euler_Angles
(
    IOobject
    (
        "Euler_Angles",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
        mesh,
dimensionedVector("Euler_Angles",dimensionSet(0, 0, 0, 0, 0),vector::zero)//,
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

volTensorField dsigmadt(
    IOobject(
        "dsigmadt",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedTensor("dsigmadt", dimensionSet(1, -1, -3, 0, 0), tensor::zero)//,
    // zeroGradientFvPatchTensorField::typeName
    );

    // volTensorField dsigmadtDIFF(
    // IOobject(
    //     "dsigmadtDIFF",
    //     runTime.timeName(),
    //     mesh,
    //     IOobject::NO_READ,
    //     IOobject::AUTO_WRITE),
    // mesh,
    // dimensionedTensor("dsigmadtDIFF", dimensionSet(1, -1, -3, 0, 0), tensor::zero)//,
    // // zeroGradientFvPatchTensorField::typeName
    // );

// volTensorField dsigmadt_crystal
// (
//     IOobject
//     (
//         "dsigmadt_crystal",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//         mesh,
// dimensionedTensor("dsigmadt_crystal",dimensionSet(1, -1, -3, 0, 0),tensor::zero),
//                 zeroGradientFvPatchTensorField::typeName
// );

#include "createPhi.H"



// volTensorField test_rot_diff
// (
//     IOobject
//     (
//         "test_rot_diff",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::AUTO_WRITE
//     ),
//         mesh,
// dimensionedTensor("test_rot_diff",dimensionSet(0, 0, 0, 0, 0),(tensor::zero))//,
//                 //zeroGradientFvPatchTensorField::typeName
// );

// PtrList<volTensorField> test_rot_diff(N_Ori); // pointer to phase fields

// for (label count = 0; count < N_Ori; count++) // how to assign the popbal array of fields
// {

//     word name_ni("testRot." + name(count));

//     test_rot_diff.set(
//         count,
//         new volTensorField(
//             IOobject(
//                 name_ni,
//                 runTime.timeName(),
//                 mesh,
//                 IOobject::READ_IF_PRESENT, // IOobject::MUST_READ,
//                 IOobject::AUTO_WRITE),
//             mesh,
//             dimensionedTensor("ni", dimensionSet(0, 0, 0, 0, 0), (tensor::zero)),
//             zeroGradientFvPatchTensorField::typeName));
// }


//   forAll(mesh.C(), cell)
//   {

//     for (label h = 0; h < N_Ori; h++)
//     {

//       double phi_1_total = phi1[h]; // ROT_1[cell];// phi1[h] ;//+ omega_tot[cell].component(5);

//       double phi_2_total = phi2[h]; // ROT_2[cell];// phi2[h] ;//+ omega_tot[cell].component(2);

//       double phi_3_total = phi3[h];

//       Fastor::Tensor<double, 3, 3> Rot_matrix(0.0);

//       Rot_matrix(0, 0) = Foam::cos(phi_2_total) * Foam::cos(phi_3_total);
//       Rot_matrix(0, 1) = (Foam::sin(phi_1_total) * Foam::sin(phi_2_total) * Foam::cos(phi_3_total)) - (Foam::cos(phi_1_total) * Foam::sin(phi_3_total));
//       Rot_matrix(0, 2) = (Foam::cos(phi_1_total) * Foam::sin(phi_2_total) * Foam::cos(phi_3_total)) + (Foam::sin(phi_1_total) * Foam::sin(phi_3_total));

//       Rot_matrix(1, 0) = Foam::cos(phi_2_total) * Foam::sin(phi_3_total);
//       Rot_matrix(1, 1) = (Foam::sin(phi_1_total) * Foam::sin(phi_2_total) * Foam::sin(phi_3_total)) + (Foam::cos(phi_1_total) * Foam::cos(phi_3_total));
//       Rot_matrix(1, 2) = (Foam::cos(phi_1_total) * Foam::sin(phi_2_total) * Foam::sin(phi_3_total)) - (Foam::sin(phi_1_total) * Foam::cos(phi_3_total));

//       Rot_matrix(2, 0) = -Foam::sin(phi_2_total);
//       Rot_matrix(2, 1) = Foam::sin(phi_1_total) * Foam::cos(phi_2_total);
//       Rot_matrix(2, 2) = Foam::cos(phi_1_total) * Foam::cos(phi_2_total);


//     // if(cell==240){
//     //     Info<<h<<endl;
//     // Fastor::print("Rot: ",Rot_matrix);
//     // }


//       test_rot_diff[h][cell].component(0)=Rot_matrix(0, 0);
//       test_rot_diff[h][cell].component(1)=Rot_matrix(0, 1);
//       test_rot_diff[h][cell].component(2)=Rot_matrix(0, 2);
//       test_rot_diff[h][cell].component(3)=Rot_matrix(1, 0);
//       test_rot_diff[h][cell].component(4)=Rot_matrix(1, 1);
//       test_rot_diff[h][cell].component(5)=Rot_matrix(1, 2);
//       test_rot_diff[h][cell].component(6)=Rot_matrix(2, 0);
//       test_rot_diff[h][cell].component(7)=Rot_matrix(2, 1);
//       test_rot_diff[h][cell].component(8)=Rot_matrix(2, 2);

//     }

//   }



// Fastor::Tensor<double, 3> TEST1(0.0);
// Fastor::Tensor<double, 3> TEST2(0.0);

// TEST1(0) = 2.0;
// TEST1(1) = 4.0;
// TEST1(2) = 5.0;

// TEST2(0) = 3.0;
// TEST2(1) = 7.0;
// TEST2(2) = 8.0;


// Fastor::print(TEST1);
// Fastor::print(TEST2);
// Fastor::print(outer(TEST1,TEST2));


// exit(0);


// forAll(mesh.C(), cell)
//   {
//     // test_rot_diff[cell].component(0)=1.0;
//     // test_rot_diff[cell].component(4)=2.0;
//     // test_rot_diff[cell].component(8)=3.0;

//     test_rot_diff[0][cell].component(0)=1.0/Foam::sqrt(2.0);
//     test_rot_diff[0][cell].component(1)=1.0/Foam::sqrt(6.0);
//     test_rot_diff[0][cell].component(2)=1.0/Foam::sqrt(3.0);
//     test_rot_diff[0][cell].component(3)=-1.0/Foam::sqrt(2.0);
//     test_rot_diff[0][cell].component(4)=1.0/Foam::sqrt(6.0);
//     test_rot_diff[0][cell].component(5)=1.0/Foam::sqrt(3.0);
//     test_rot_diff[0][cell].component(6)=0.0;
//     test_rot_diff[0][cell].component(7)=-2.0/Foam::sqrt(6.0);
//     test_rot_diff[0][cell].component(8)=1.0/Foam::sqrt(3.0);
//   }



// Fastor::Tensor<double, 3, 3> TEST_1(0.0);

// TEST_1(0,0) = -1.0;
// TEST_1(0,1) = 2.0;
// TEST_1(0,2) = 3.0;
// TEST_1(1,0) = 4.0;
// TEST_1(1,1) = -5.0;
// TEST_1(1,2) = 6.0;
// TEST_1(2,0) = 7.0;
// TEST_1(2,1) = 8.0;
// TEST_1(2,2) = 9.0;

// Fastor::Tensor<double, 3, 3> TEST_2(0.0);

// TEST_2(0,0) = 2.0;
// TEST_2(0,1) = -3.0;
// TEST_2(0,2) = -4.0;
// TEST_2(1,0) = 5.0;
// TEST_2(1,1) = -6.0;
// TEST_2(1,2) = 7.0;
// TEST_2(2,0) = 8.0;
// TEST_2(2,1) = -9.0;
// TEST_2(2,2) = 10.0;

// Fastor::print(Fastor::matmul(TEST_1,TEST_2));



// tensor T1(-1.0,2.0,3.0,4.0,-5.0,6.0,7.0,8.0,9.0);
// tensor T2(2.0,-3.0,-4.0,5.0,-6.0,7.0,8.0,-9.0,10.0);

// Info<<"Foam dot "<<(T1 & T2);

// exit(0);


// tensor T3(1.0,0.0,0.0,0.0,0.9397,-0.3420,0.0,0.3420,0.9397);
// tensor T4(1.0,0.0,0.0,0.0,2.0,0.0,0.0,0.0,3.0);

// Info<<"Foam dot "<<((T3 & T4) & T3.T());

// exit(0);


// //// TRYING TO MAKE SENSE OF ANGLES//
    // double phi_1_total1=phi1[1];//0.0314159;//0.5*M_PI;
    // double phi_2_total1=phi2[1];
    // double phi_3_total1=phi3[1];

    //   Fastor::Tensor<double, 3, 3> Rot_matrixtest1(0.0);

    //   Rot_matrixtest1(0, 0) = Foam::cos(phi_2_total1) * Foam::cos(phi_3_total1);
    //   Rot_matrixtest1(0, 1) = (Foam::sin(phi_1_total1) * Foam::sin(phi_2_total1) * Foam::cos(phi_3_total1)) - (Foam::cos(phi_1_total1) * Foam::sin(phi_3_total1));
    //   Rot_matrixtest1(0, 2) = (Foam::cos(phi_1_total1) * Foam::sin(phi_2_total1) * Foam::cos(phi_3_total1)) + (Foam::sin(phi_1_total1) * Foam::sin(phi_3_total1));

    //   Rot_matrixtest1(1, 0) = Foam::cos(phi_2_total1) * Foam::sin(phi_3_total1);
    //   Rot_matrixtest1(1, 1) = (Foam::sin(phi_1_total1) * Foam::sin(phi_2_total1) * Foam::sin(phi_3_total1)) + (Foam::cos(phi_1_total1) * Foam::cos(phi_3_total1));
    //   Rot_matrixtest1(1, 2) = (Foam::cos(phi_1_total1) * Foam::sin(phi_2_total1) * Foam::sin(phi_3_total1)) - (Foam::sin(phi_1_total1) * Foam::cos(phi_3_total1));

    //   Rot_matrixtest1(2, 0) = -Foam::sin(phi_2_total1);
    //   Rot_matrixtest1(2, 1) = Foam::sin(phi_1_total1) * Foam::cos(phi_2_total1);
    //   Rot_matrixtest1(2, 2) = Foam::cos(phi_1_total1) * Foam::cos(phi_2_total1);



    // double phi_1_total2=phi1[1]-0.7*M_PI;//0.5*M_PI;//0.0314159;//0.5*M_PI;//after 1e-7
    // double phi_2_total2=phi2[1];
    // double phi_3_total2=phi3[1];

    //   Fastor::Tensor<double, 3, 3> Rot_matrixtest(0.0);

    //   Rot_matrixtest(0, 0) = Foam::cos(phi_2_total2) * Foam::cos(phi_3_total2);
    //   Rot_matrixtest(0, 1) = (Foam::sin(phi_1_total2) * Foam::sin(phi_2_total2) * Foam::cos(phi_3_total2)) - (Foam::cos(phi_1_total2) * Foam::sin(phi_3_total2));
    //   Rot_matrixtest(0, 2) = (Foam::cos(phi_1_total2) * Foam::sin(phi_2_total2) * Foam::cos(phi_3_total2)) + (Foam::sin(phi_1_total2) * Foam::sin(phi_3_total2));

    //   Rot_matrixtest(1, 0) = Foam::cos(phi_2_total2) * Foam::sin(phi_3_total2);
    //   Rot_matrixtest(1, 1) = (Foam::sin(phi_1_total2) * Foam::sin(phi_2_total2) * Foam::sin(phi_3_total2)) + (Foam::cos(phi_1_total2) * Foam::cos(phi_3_total2));
    //   Rot_matrixtest(1, 2) = (Foam::cos(phi_1_total2) * Foam::sin(phi_2_total2) * Foam::sin(phi_3_total2)) - (Foam::sin(phi_1_total2) * Foam::cos(phi_3_total2));

    //   Rot_matrixtest(2, 0) = -Foam::sin(phi_2_total2);
    //   Rot_matrixtest(2, 1) = Foam::sin(phi_1_total2) * Foam::cos(phi_2_total2);
    //   Rot_matrixtest(2, 2) = Foam::cos(phi_1_total2) * Foam::cos(phi_2_total2);

    // Fastor::print("rot matrix initial: ",Rot_matrixtest1);
    // Fastor::print("rot matrix aftre - pi/2: ",Rot_matrixtest);

    //   exit(0);