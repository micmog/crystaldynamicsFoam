{


  // gradU = fvc::grad(U);
  // gradU.correctBoundaryConditions();

  volTensorField TgradU(fvc::grad(U)().T());
  TgradU.correctBoundaryConditions(); // already called by line above

  omega_dot = 0.5 * (TgradU - fvc::grad(U));   
omega_dot.correctBoundaryConditions();
  //NOTE!!! Openfoam defines grad in the opposite manner to solid mechanics textbooks, it's actually grad transpose that we want so this has to be transposed again or just ahev a - sign

  epsilon_dot = 0.5 * (TgradU + fvc::grad(U));//- (1.0/3.0)*(fvc::div(phi))*I; // this looks like dev(gradU)
epsilon_dot.correctBoundaryConditions();
  // dsigmadt = -(sigma & omega_dot) + (omega_dot & sigma) ;//-(tr(epsilon_dot)*sigma);;// 

dsigmadt*=0.0;



//
// gradU = TgradU;//fvc::grad(U);
// L = epsilon_dot + omega_dot;
//

GBREF*=0.0;
  for (label h = 0; h < N_Ori; h++)
  {
    Omega_P[h] *= 0.0;
    volScalarField & ni = PopBal[h];

    GBREF+=ni*ni;

    PopBal_av[h] = fvc::average(ni); //takes the average from all cells sharing a face - this way the information about phases in adjacent cells is shared with the local cell
    PopBal_av[h].correctBoundaryConditions();
  }
  GBREF.correctBoundaryConditions();

  D_P *= 0.0;

  p *= 0.0;

  














  for (label i = 0; i < N_Ori; i++)
  {
    for (label countS = 0; countS < MatList[Phase2Material[i]].N_Slips; countS++)
    {
      // g_alpha_dots[i][countS].write();
      gamma_alpha_dots[i][countS]*=0.0;
    }
  }





// Euler_Angles*=0.0;

  forAll(mesh.C(), cell)
  {



    // label active_phases[N_Ori] =
    // {
    //     0
    // };
    // label n_active_phases = 0;



    // for (label count = 0; count < N_Ori; count++)
    // {
        
    //     if (PopBal_av[count][cell] > 0.0)
    //     {
    //         active_phases[n_active_phases++] = count;
    //     }

    // }




    Fastor::Tensor<double, 3, 3> Fsigma(0.0); // CAST OPENFOAM VARIABLES INTO FASTOR VARIABLES
    Fastor::Tensor<double, 3, 3> Fepsilon_dot(0.0);

    Fastor::Tensor<double, 3, 3> F_D_P(0.0);//plastic strain rate
    Fastor::Tensor<double, 3, 3> F_Omega_P(0.0);//plastic rotation rate

    int count_index = 0;
    for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++)
      {
        // Info<<i<<","<<j<<endl;
        Fepsilon_dot(i, j) = epsilon_dot[cell].component(count_index);
        Fsigma(i, j) = sigtemp[cell].component(count_index);


        count_index++;
      }
    }

    Fastor::Tensor<double, 3, 3, 3, 3> CprimeSum(0.0); // global stiffness at this cell

          Fastor::Tensor<double, 3, 3> P_alpha(0.0);
          Fastor::Tensor<double, 3, 3> omega_alpha(0.0);


      Fastor::Tensor<double, 3, 3> Rot_matrix(0.0);

      Rot_matrix(0, 0) = RotationMatrix[cell].component(0);//Foam::cos(phi_2_total) * Foam::cos(phi_3_total);
      Rot_matrix(0, 1) = RotationMatrix[cell].component(1);//(Foam::sin(phi_1_total) * Foam::sin(phi_2_total) * Foam::cos(phi_3_total)) - (Foam::cos(phi_1_total) * Foam::sin(phi_3_total));
      Rot_matrix(0, 2) = RotationMatrix[cell].component(2);//(Foam::cos(phi_1_total) * Foam::sin(phi_2_total) * Foam::cos(phi_3_total)) + (Foam::sin(phi_1_total) * Foam::sin(phi_3_total));

      Rot_matrix(1, 0) = RotationMatrix[cell].component(3);//Foam::cos(phi_2_total) * Foam::sin(phi_3_total);
      Rot_matrix(1, 1) = RotationMatrix[cell].component(4);//(Foam::sin(phi_1_total) * Foam::sin(phi_2_total) * Foam::sin(phi_3_total)) + (Foam::cos(phi_1_total) * Foam::cos(phi_3_total));
      Rot_matrix(1, 2) = RotationMatrix[cell].component(5);//(Foam::cos(phi_1_total) * Foam::sin(phi_2_total) * Foam::sin(phi_3_total)) - (Foam::sin(phi_1_total) * Foam::cos(phi_3_total));

      Rot_matrix(2, 0) = RotationMatrix[cell].component(6);//-Foam::sin(phi_2_total);
      Rot_matrix(2, 1) = RotationMatrix[cell].component(7);//Foam::sin(phi_1_total) * Foam::cos(phi_2_total);
      Rot_matrix(2, 2) = RotationMatrix[cell].component(8);//Foam::cos(phi_1_total) * Foam::cos(phi_2_total);



    for (label h = 0; h < N_Ori; h++)
    {

      Fastor::Tensor<double, 3, 3, 3, 3> Cprime(0.0); // 


      // if (PopBal[h][cell] > SMALL)
      // { // if the phase is there, work out the stiffness contribution

        for (int i = 0; i < 3; i++)
        {
          for (int j = 0; j < 3; j++)
          {
            for (int k = 0; k < 3; k++)
            {
              for (int l = 0; l < 3; l++)
              {

                double value = 0.0;
                for (int m = 0; m < 3; m++)
                {
                  for (int n = 0; n < 3; n++)
                  {
                    for (int p = 0; p < 3; p++)
                    {
                      for (int q = 0; q < 3; q++)
                      {
                        // if (Cprime2[Phase2Material[h]](m, n, p, q) > SMALL)
                        // {
                          value += Rot_matrix(i, m) * Rot_matrix(j, n) * Rot_matrix(k, p) * Rot_matrix(l, q) * Cprime2[Phase2Material[h]](m, n, p, q); // C(m,n,p,q);
                        // }
                      }
                    }
                  }
                }

                // value*=PopBal[h][cell];

                Cprime(i, j, k, l) = value;
                CprimeSum(i, j, k, l)+=PopBal[h][cell]*Cprime(i, j, k, l);
              }
            }
          }
        }

        

      // }//

        // If grain is there
        // Have the stiffness tensor now
        // need Dp and omegaP


        for (int s = 0; s < MatList[Phase2Material[h]].N_Slips; s++)
        {

          Fastor::Tensor<double, 3> F_n_alpha(0.0);
          Fastor::Tensor<double, 3> F_m_alpha(0.0);
          // Fastor::Tensor<double, 3, 3> P_alpha(0.0);
          // Fastor::Tensor<double, 3, 3> omega_alpha(0.0);

          F_n_alpha(0) = MatList[Phase2Material[h]].SP_norm[s].component(0); // n_alpha[Phase2Material[h]][s].component(0);
          F_n_alpha(1) = MatList[Phase2Material[h]].SP_norm[s].component(1); // n_alpha[Phase2Material[h]][s].component(1);
          F_n_alpha(2) = MatList[Phase2Material[h]].SP_norm[s].component(2); // n_alpha[Phase2Material[h]][s].component(2);

          F_m_alpha(0) = MatList[Phase2Material[h]].SP_direction[s].component(0); // m_alpha[Phase2Material[h]][s].component(0);
          F_m_alpha(1) = MatList[Phase2Material[h]].SP_direction[s].component(1); // m_alpha[Phase2Material[h]][s].component(1);
          F_m_alpha(2) = MatList[Phase2Material[h]].SP_direction[s].component(2); // m_alpha[Phase2Material[h]][s].component(2);


          Fastor::Tensor<double, 3> F_n_alpha_Rot(0.0);
          Fastor::Tensor<double, 3> F_m_alpha_Rot(0.0);

          for (int i = 0; i < 3; i++)
          {
            for (int j = 0; j < 3; j++)
            { // Rotated slip information
              F_n_alpha_Rot(i) += Rot_matrix(i, j) * F_n_alpha(j);
              F_m_alpha_Rot(i) += Rot_matrix(i, j) * F_m_alpha(j);
            }
          }


          P_alpha = 0.5 * (outer(F_m_alpha_Rot, F_n_alpha_Rot) + outer(F_n_alpha_Rot, F_m_alpha_Rot));
          omega_alpha = 0.5 * (outer(F_m_alpha_Rot, F_n_alpha_Rot) - outer(F_n_alpha_Rot, F_m_alpha_Rot));

          // Fastor::Tensor<double, 3, 3> Schmid_T(0.0);

          // Schmid_T = outer(F_m_alpha_Rot, F_n_alpha_Rot);
          
          double tau_alpha = 0.0;
          for (int i = 0; i < 3; i++)
          {
            for (int j = 0; j < 3; j++)
            {
              tau_alpha += Fsigma(i, j)* P_alpha(i, j);
            }
          }

          double k_exponent = MatList[Phase2Material[h]].k_exponent;         // 10.0;//12.0
          double ref_shear_rate = MatList[Phase2Material[h]].ref_shear_rate; // 1e-3;
          

          gamma_alpha_dots[h][s][cell] = ref_shear_rate * Foam::pow((mag(tau_alpha / g_alphas[h][s][cell])), k_exponent) * sign(tau_alpha); // source term for RHS of g_alpha equations



          F_D_P += PopBal[h][cell]*(P_alpha * gamma_alpha_dots[h][s][cell]);
          F_Omega_P += PopBal[h][cell]*(omega_alpha * gamma_alpha_dots[h][s][cell]);


        }         // end of loop over slip systems

        // sumgamma[cell]+=PopBal[h][cell]*SUMGAMMA;





    } // end loop over phase fields



        int OP_index = 0;
        for (int i = 0; i < 3; i++)
        {
          for (int j = 0; j < 3; j++)
          {
            Omega_P[cell].component(OP_index) = F_Omega_P(i, j);
            D_P[cell].component(OP_index) =  F_D_P(i,j);
            OP_index++;
          }
        }


      enum
        {
          I,J,K,L
        };
        Fastor::Tensor<double, 3, 3> sig = Fastor::einsum<Fastor::Index<I, J, K, L>, Fastor::Index<K, L> /*,Fastor::OIndex<I,J>*/>(CprimeSum, Fepsilon_dot - F_D_P); ////Double contraction


        dsigmadt[cell] -= ((sigtemp[cell] & omega_dot[cell]) + (omega_dot[cell] & sigtemp[cell])) ;//-(tr(epsilon_dot)*sigma);;// 
        dsigmadt[cell] += ((sigtemp[cell] & Omega_P[cell]) - (Omega_P[cell] & sigtemp[cell]));




    count_index = 0;
    for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++)
      {

        dsigmadt[cell].component(count_index) += sig(i, j); 
        count_index++;
      }
    }


  }





//////// Over Boundaries
forAll(dsigmadt.boundaryField(), patchI)
{
    if
    (
        !dsigmadt.boundaryField()[patchI].coupled()
     && dsigmadt.boundaryField()[patchI].size() > 0
    )
    {
        forAll(dsigmadt.boundaryField()[patchI], faceI)
        {





    Fastor::Tensor<double, 3, 3> Fsigma(0.0); // CAST OPENFOAM VARIABLES INTO FASTOR VARIABLES
    Fastor::Tensor<double, 3, 3> Fepsilon_dot(0.0);

    Fastor::Tensor<double, 3, 3> F_D_P(0.0);//plastic strain rate
    Fastor::Tensor<double, 3, 3> F_Omega_P(0.0);//plastic rotation rate

    int count_index = 0;
    for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++)
      {
        // Info<<i<<","<<j<<endl;
        Fepsilon_dot(i, j) = epsilon_dot.boundaryField()[patchI][faceI].component(count_index);
        Fsigma(i, j) = sigtemp.boundaryField()[patchI][faceI].component(count_index);


        count_index++;
      }
    }

    Fastor::Tensor<double, 3, 3, 3, 3> CprimeSum(0.0); // global stiffness at this cell

          Fastor::Tensor<double, 3, 3> P_alpha(0.0);
          Fastor::Tensor<double, 3, 3> omega_alpha(0.0);


      Fastor::Tensor<double, 3, 3> Rot_matrix(0.0);

      Rot_matrix(0, 0) = RotationMatrix.boundaryField()[patchI][faceI].component(0);//Foam::cos(phi_2_total) * Foam::cos(phi_3_total);
      Rot_matrix(0, 1) = RotationMatrix.boundaryField()[patchI][faceI].component(1);//(Foam::sin(phi_1_total) * Foam::sin(phi_2_total) * Foam::cos(phi_3_total)) - (Foam::cos(phi_1_total) * Foam::sin(phi_3_total));
      Rot_matrix(0, 2) = RotationMatrix.boundaryField()[patchI][faceI].component(2);//(Foam::cos(phi_1_total) * Foam::sin(phi_2_total) * Foam::cos(phi_3_total)) + (Foam::sin(phi_1_total) * Foam::sin(phi_3_total));

      Rot_matrix(1, 0) = RotationMatrix.boundaryField()[patchI][faceI].component(3);//Foam::cos(phi_2_total) * Foam::sin(phi_3_total);
      Rot_matrix(1, 1) = RotationMatrix.boundaryField()[patchI][faceI].component(4);//(Foam::sin(phi_1_total) * Foam::sin(phi_2_total) * Foam::sin(phi_3_total)) + (Foam::cos(phi_1_total) * Foam::cos(phi_3_total));
      Rot_matrix(1, 2) = RotationMatrix.boundaryField()[patchI][faceI].component(5);//(Foam::cos(phi_1_total) * Foam::sin(phi_2_total) * Foam::sin(phi_3_total)) - (Foam::sin(phi_1_total) * Foam::cos(phi_3_total));

      Rot_matrix(2, 0) = RotationMatrix.boundaryField()[patchI][faceI].component(6);//-Foam::sin(phi_2_total);
      Rot_matrix(2, 1) = RotationMatrix.boundaryField()[patchI][faceI].component(7);//Foam::sin(phi_1_total) * Foam::cos(phi_2_total);
      Rot_matrix(2, 2) = RotationMatrix.boundaryField()[patchI][faceI].component(8);//Foam::cos(phi_1_total) * Foam::cos(phi_2_total);



    for (label h = 0; h < N_Ori; h++)
    {

      Fastor::Tensor<double, 3, 3, 3, 3> Cprime(0.0); // 


      // if (PopBal[h][cell] > SMALL)
      // { // if the phase is there, work out the stiffness contribution

        for (int i = 0; i < 3; i++)
        {
          for (int j = 0; j < 3; j++)
          {
            for (int k = 0; k < 3; k++)
            {
              for (int l = 0; l < 3; l++)
              {

                double value = 0.0;
                for (int m = 0; m < 3; m++)
                {
                  for (int n = 0; n < 3; n++)
                  {
                    for (int p = 0; p < 3; p++)
                    {
                      for (int q = 0; q < 3; q++)
                      {
                        // if (Cprime2[Phase2Material[h]](m, n, p, q) > SMALL)
                        // {
                          value += Rot_matrix(i, m) * Rot_matrix(j, n) * Rot_matrix(k, p) * Rot_matrix(l, q) * Cprime2[Phase2Material[h]](m, n, p, q); // C(m,n,p,q);
                        // }
                      }
                    }
                  }
                }

                // value*=PopBal[h][cell];

                Cprime(i, j, k, l) = value;
                CprimeSum(i, j, k, l)+=PopBal[h].boundaryField()[patchI][faceI]*Cprime(i, j, k, l);
              }
            }
          }
        }

        

    //   // }//

    //     // If grain is there
    //     // Have the stiffness tensor now
    //     // need Dp and omegaP


        for (int s = 0; s < MatList[Phase2Material[h]].N_Slips; s++)
        {

          Fastor::Tensor<double, 3> F_n_alpha(0.0);
          Fastor::Tensor<double, 3> F_m_alpha(0.0);

          F_n_alpha(0) = MatList[Phase2Material[h]].SP_norm[s].component(0); // n_alpha[Phase2Material[h]][s].component(0);
          F_n_alpha(1) = MatList[Phase2Material[h]].SP_norm[s].component(1); // n_alpha[Phase2Material[h]][s].component(1);
          F_n_alpha(2) = MatList[Phase2Material[h]].SP_norm[s].component(2); // n_alpha[Phase2Material[h]][s].component(2);

          F_m_alpha(0) = MatList[Phase2Material[h]].SP_direction[s].component(0); // m_alpha[Phase2Material[h]][s].component(0);
          F_m_alpha(1) = MatList[Phase2Material[h]].SP_direction[s].component(1); // m_alpha[Phase2Material[h]][s].component(1);
          F_m_alpha(2) = MatList[Phase2Material[h]].SP_direction[s].component(2); // m_alpha[Phase2Material[h]][s].component(2);


          Fastor::Tensor<double, 3> F_n_alpha_Rot(0.0);
          Fastor::Tensor<double, 3> F_m_alpha_Rot(0.0);

          for (int i = 0; i < 3; i++)
          {
            for (int j = 0; j < 3; j++)
            { // Rotated slip information
              F_n_alpha_Rot(i) += Rot_matrix(i, j) * F_n_alpha(j);
              F_m_alpha_Rot(i) += Rot_matrix(i, j) * F_m_alpha(j);
            }
          }


          P_alpha = 0.5 * (outer(F_m_alpha_Rot, F_n_alpha_Rot) + outer(F_n_alpha_Rot, F_m_alpha_Rot));
          omega_alpha = 0.5 * (outer(F_m_alpha_Rot, F_n_alpha_Rot) - outer(F_n_alpha_Rot, F_m_alpha_Rot));

          // Fastor::Tensor<double, 3, 3> Schmid_T(0.0);

          // Schmid_T = outer(F_m_alpha_Rot, F_n_alpha_Rot);
          
          double tau_alpha = 0.0;
          for (int i = 0; i < 3; i++)
          {
            for (int j = 0; j < 3; j++)
            {
              tau_alpha += Fsigma(i, j)* P_alpha(i, j);
            }
          }

          double k_exponent = MatList[Phase2Material[h]].k_exponent;         // 10.0;//12.0
          double ref_shear_rate = MatList[Phase2Material[h]].ref_shear_rate; // 1e-3;
          

          gamma_alpha_dots[h][s].boundaryFieldRef()[patchI][faceI] = ref_shear_rate * Foam::pow((mag(tau_alpha / g_alphas[h][s].boundaryField()[patchI][faceI])), k_exponent) * sign(tau_alpha); // source term for RHS of g_alpha equations



          F_D_P += PopBal[h].boundaryField()[patchI][faceI]*(P_alpha * gamma_alpha_dots[h][s].boundaryField()[patchI][faceI]);
          F_Omega_P += PopBal[h].boundaryField()[patchI][faceI]*(omega_alpha * gamma_alpha_dots[h][s].boundaryField()[patchI][faceI]);


        }         // end of loop over slip systems

    //     // sumgamma[cell]+=PopBal[h][cell]*SUMGAMMA;





    } // end loop over phase fields



        int OP_index = 0;
        for (int i = 0; i < 3; i++)
        {
          for (int j = 0; j < 3; j++)
          {
            Omega_P.boundaryFieldRef()[patchI][faceI].component(OP_index) = F_Omega_P(i, j);
            D_P.boundaryFieldRef()[patchI][faceI].component(count_index) =  F_D_P(i,j);
            OP_index++;
          }
        }


      enum
        {
          I,J,K,L
        };
        Fastor::Tensor<double, 3, 3> sig = Fastor::einsum<Fastor::Index<I, J, K, L>, Fastor::Index<K, L> /*,Fastor::OIndex<I,J>*/>(CprimeSum, Fepsilon_dot - F_D_P); ////Double contraction


        dsigmadt.boundaryFieldRef()[patchI][faceI] -= ((sigtemp.boundaryField()[patchI][faceI] & omega_dot.boundaryField()[patchI][faceI]) + (omega_dot.boundaryField()[patchI][faceI] & sigtemp.boundaryField()[patchI][faceI])) ;//-(tr(epsilon_dot)*sigma);;// 
        dsigmadt.boundaryFieldRef()[patchI][faceI] += ((sigtemp.boundaryField()[patchI][faceI] & Omega_P.boundaryField()[patchI][faceI]) - (Omega_P.boundaryField()[patchI][faceI] & sigtemp.boundaryField()[patchI][faceI]));




    count_index = 0;
    for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++)
      {

        dsigmadt.boundaryFieldRef()[patchI][faceI].component(count_index) += sig(i, j); 
        count_index++;
      }
    }



        
        }
    }
}







// // Over Boundaries
// forAll(dsigmadt.boundaryField(), patchI)
// {
//     if
//     (
//         !dsigmadt.boundaryField()[patchI].coupled()
//      && dsigmadt.boundaryField()[patchI].size() > 0
//     )
//     {
//         forAll(dsigmadt.boundaryField()[patchI], faceI)
//         {


//     Fastor::Tensor<double, 3, 3> Fsigma(0.0); // CAST OPENFOAM VARIABLES INTO FASTOR VARIABLES
//     Fastor::Tensor<double, 3, 3> Fepsilon_dot(0.0);
//     // Fastor::Tensor<double, 3, 3> Fsigmaav(0.0);
//     int count_index = 0;
//     for (int i = 0; i < 3; i++)
//     {
//       for (int j = 0; j < 3; j++)
//       {
//         // Info<<i<<","<<j<<endl;
//         Fepsilon_dot(i, j) = epsilon_dot.boundaryField()[patchI][faceI].component(count_index);
//         Fsigma(i, j) = sigtemp.boundaryField()[patchI][faceI].component(count_index);
//         // Fsigmaav(i, j) = sigmaav[cell].component(count_index);

//         count_index++;
//       }
//     }


//     Fastor::Tensor<double, 3, 3> sigAv(0.0);
//     Fastor::Tensor<double, 3, 3> D_P_Average(0.0);


//         label largest_field=-1;
//     double current_largest_value = -1.0;

//     sumgamma.boundaryFieldRef()[patchI][faceI]*=0.0;
//     for (label h = 0; h < N_Ori; h++)
//     {


//         if(PopBal[h].boundaryField()[patchI][faceI]>current_largest_value){
//         current_largest_value = PopBal[h].boundaryField()[patchI][faceI];
      
//         largest_field = h;
//       }



//        Fastor::Tensor<double, 3, 3, 3, 3> Cprime(0.0); // global stiffness at this cell

//       double phi_1_total = Rot1[h].boundaryField()[patchI][faceI]; // ROT_1[cell];// phi1[h] ;//+ omega_tot[cell].component(5);

//       double phi_2_total = Rot2[h].boundaryField()[patchI][faceI]; // ROT_2[cell];// phi2[h] ;//+ omega_tot[cell].component(2);

//       double phi_3_total = Rot3[h].boundaryField()[patchI][faceI]; // ROT_3[cell];// phi3[h] ;//+ omega_tot[cell].component(1);

//           // if(PopBal[h].boundaryField()[patchI][faceI]>0.5){

//     // }

//     Fastor::Tensor<double, 3, 3> Rot_matrix(0.0);

//       Rot_matrix(0, 0) = Foam::cos(phi_2_total) * Foam::cos(phi_3_total);
//       Rot_matrix(0, 1) = (Foam::sin(phi_1_total) * Foam::sin(phi_2_total) * Foam::cos(phi_3_total)) - (Foam::cos(phi_1_total) * Foam::sin(phi_3_total));
//       Rot_matrix(0, 2) = (Foam::cos(phi_1_total) * Foam::sin(phi_2_total) * Foam::cos(phi_3_total)) + (Foam::sin(phi_1_total) * Foam::sin(phi_3_total));

//       Rot_matrix(1, 0) = Foam::cos(phi_2_total) * Foam::sin(phi_3_total);
//       Rot_matrix(1, 1) = (Foam::sin(phi_1_total) * Foam::sin(phi_2_total) * Foam::sin(phi_3_total)) + (Foam::cos(phi_1_total) * Foam::cos(phi_3_total));
//       Rot_matrix(1, 2) = (Foam::cos(phi_1_total) * Foam::sin(phi_2_total) * Foam::sin(phi_3_total)) - (Foam::sin(phi_1_total) * Foam::cos(phi_3_total));

//       Rot_matrix(2, 0) = -Foam::sin(phi_2_total);
//       Rot_matrix(2, 1) = Foam::sin(phi_1_total) * Foam::cos(phi_2_total);
//       Rot_matrix(2, 2) = Foam::cos(phi_1_total) * Foam::cos(phi_2_total);
//       // https://en.wikipedia.org/wiki/Rotation_matrix -
//       //"represents an extrinsic rotation whose (improper) Euler angles are α, β, γ, about axes x, y, z" - phi1 ==alpha, phi2==beta,phi3==gamma



//       if (PopBal[h].boundaryField()[patchI][faceI] > SMALL)
//       { // if the phase is there, work out the stiffness contribution

//         for (int i = 0; i < 3; i++)
//         {
//           for (int j = 0; j < 3; j++)
//           {
//             for (int k = 0; k < 3; k++)
//             {
//               for (int l = 0; l < 3; l++)
//               {

//                 double value = 0.0;
//                 for (int m = 0; m < 3; m++)
//                 {
//                   for (int n = 0; n < 3; n++)
//                   {
//                     for (int p = 0; p < 3; p++)
//                     {
//                       for (int q = 0; q < 3; q++)
//                       {
//                         // if (Cprime2[Phase2Material[h]](m, n, p, q) > SMALL)
//                         // {
//                           value += Rot_matrix(i, m) * Rot_matrix(j, n) * Rot_matrix(k, p) * Rot_matrix(l, q) * Cprime2[Phase2Material[h]](m, n, p, q); // C(m,n,p,q);
//                         // }
//                       }
//                     }
//                   }
//                 }

//                 // value*=PopBal[h][cell];

//                 Cprime(i, j, k, l) = value;
//               }
//             }
//           }
//         }

//       }//


// // If grain is there
//         // Have the stiffness tensor now
//         // need Dp and omegaP

//         Fastor::Tensor<double, 3, 3> F_D_P(0.0);
//         Fastor::Tensor<double, 3, 3> F_Omega_P(0.0);

//       double SUMGAMMA=0.0;
//       for (int s = 0; s < MatList[Phase2Material[h]].N_Slips; s++)
//         {

//           Fastor::Tensor<double, 3> F_n_alpha(0.0);
//           Fastor::Tensor<double, 3> F_m_alpha(0.0);
//           Fastor::Tensor<double, 3, 3> P_alpha(0.0);
//           Fastor::Tensor<double, 3, 3> omega_alpha(0.0);

//           F_n_alpha(0) = MatList[Phase2Material[h]].SP_norm[s].component(0); // n_alpha[Phase2Material[h]][s].component(0);
//           F_n_alpha(1) = MatList[Phase2Material[h]].SP_norm[s].component(1); // n_alpha[Phase2Material[h]][s].component(1);
//           F_n_alpha(2) = MatList[Phase2Material[h]].SP_norm[s].component(2); // n_alpha[Phase2Material[h]][s].component(2);

//           F_m_alpha(0) = MatList[Phase2Material[h]].SP_direction[s].component(0); // m_alpha[Phase2Material[h]][s].component(0);
//           F_m_alpha(1) = MatList[Phase2Material[h]].SP_direction[s].component(1); // m_alpha[Phase2Material[h]][s].component(1);
//           F_m_alpha(2) = MatList[Phase2Material[h]].SP_direction[s].component(2); // m_alpha[Phase2Material[h]][s].component(2);


//           Fastor::Tensor<double, 3> F_n_alpha_Rot(0.0);
//           Fastor::Tensor<double, 3> F_m_alpha_Rot(0.0);

//           for (int i = 0; i < 3; i++)
//           {
//             for (int j = 0; j < 3; j++)
//             { // Rotated slip information
//               F_n_alpha_Rot(i) += Rot_matrix(i, j) * F_n_alpha(j);
//               F_m_alpha_Rot(i) += Rot_matrix(i, j) * F_m_alpha(j);
//             }
//           }


//           P_alpha = 0.5 * (outer(F_m_alpha_Rot, F_n_alpha_Rot) + outer(F_n_alpha_Rot, F_m_alpha_Rot));
//           omega_alpha = 0.5 * (outer(F_m_alpha_Rot, F_n_alpha_Rot) - outer(F_n_alpha_Rot, F_m_alpha_Rot));

//           Fastor::Tensor<double, 3, 3> Schmid_T(0.0);

//           Schmid_T = outer(F_m_alpha_Rot, F_n_alpha_Rot);
//           double tau_alpha = 0.0;
//           for (int i = 0; i < 3; i++)
//           {
//             for (int j = 0; j < 3; j++)
//             {
//               // tau_alpha += Fsigma(i, j)* Schmid_T(i, j);//Fsigmaav(i,j) 
//               tau_alpha += Fsigma(i, j)* P_alpha(i, j);//Fsigmaav(i,j) 
//             }
//           }

//           double k_exponent = MatList[Phase2Material[h]].k_exponent;         // 10.0;//12.0
//           double ref_shear_rate = MatList[Phase2Material[h]].ref_shear_rate; // 1e-3;

//           gamma_alpha_dots[h][s].boundaryFieldRef()[patchI][faceI] = /*Foam::pow(PopBal[h].boundaryField()[patchI][faceI],10.0) */ ref_shear_rate * Foam::pow((mag(tau_alpha / g_alphas[h][s].boundaryField()[patchI][faceI])), k_exponent) * sign(tau_alpha); // source term for RHS of g_alpha equations


//           // sumgamma.boundaryFieldRef()[patchI][faceI]+= mag(gamma_alpha_dots[h][s].boundaryField()[patchI][faceI]);//for next time-step

//           SUMGAMMA+= mag(gamma_alpha_dots[h][s].boundaryField()[patchI][faceI]);//for next time-step

//           F_D_P += (P_alpha * gamma_alpha_dots[h][s].boundaryField()[patchI][faceI]);
//           F_Omega_P += (omega_alpha * gamma_alpha_dots[h][s].boundaryField()[patchI][faceI]);
//           // }

//         } // end of loop over slip systems
//         sumgamma.boundaryFieldRef()[patchI][faceI]+=PopBal[h].boundaryField()[patchI][faceI]*SUMGAMMA;

//         int OP_index = 0;
//         for (int i = 0; i < 3; i++)
//         {
//           for (int j = 0; j < 3; j++)
//           {
//             Omega_P[h].boundaryFieldRef()[patchI][faceI].component(OP_index) = F_Omega_P(i, j); //Philip, how to assign the boundary face values
//             OP_index++;
//           }
//         }

//         enum
//         {
//           I,J,K,L
//         };
//         Fastor::Tensor<double, 3, 3> sig = Fastor::einsum<Fastor::Index<I, J, K, L>, Fastor::Index<K, L> /*,Fastor::OIndex<I,J>*/>(Cprime, Fepsilon_dot - F_D_P); ////Double contraction
//         // Fastor::Tensor<double, 3, 3> sig = 80e9*(Fepsilon_dot - F_D_P);

//         sigAv += (PopBal[h].boundaryField()[patchI][faceI]) * sig;///Sum_Phases[cell];
//         // Omega_P_Average+=PopBal[h][cell]*F_Omega_P;
//       // }
//         D_P_Average += PopBal[h].boundaryField()[patchI][faceI]  * F_D_P;

//         Fastor::Tensor<double, 3, 3> F_D_E(0.0);

//         F_D_E = Fepsilon_dot-F_D_P;

// ///////Sum_Phases.boundaryField()[patchI][faceI]

//           dsigmadt.boundaryFieldRef()[patchI][faceI] -=(PopBal[h].boundaryField()[patchI][faceI]) *((sigtemp.boundaryField()[patchI][faceI] & omega_dot.boundaryField()[patchI][faceI]) + (omega_dot.boundaryField()[patchI][faceI] & sigtemp.boundaryField()[patchI][faceI]));
//           dsigmadt.boundaryFieldRef()[patchI][faceI] += (PopBal[h].boundaryField()[patchI][faceI]) *((sigtemp.boundaryField()[patchI][faceI] & Omega_P[h].boundaryField()[patchI][faceI]) - (Omega_P[h].boundaryField()[patchI][faceI] & sigtemp.boundaryField()[patchI][faceI]));//Philip, can we assign boundary fields with dot products in?
//           // dsigmadt.boundaryFieldRef()[patchI][faceI] -= (PopBal[h].boundaryField()[patchI][faceI]) *(sigtemp.boundaryField()[patchI][faceI]*(F_D_E(0,0)+F_D_E(1,1)+F_D_E(2,2)));



//     }// end loop over phase fields

//       // Euler_Angles.boundaryFieldRef()[patchI][faceI].component(0) = Rot1[largest_field].boundaryField()[patchI][faceI];
//       // Euler_Angles.boundaryFieldRef()[patchI][faceI].component(1) = Rot2[largest_field].boundaryField()[patchI][faceI];
//       // Euler_Angles.boundaryFieldRef()[patchI][faceI].component(2) = Rot3[largest_field].boundaryField()[patchI][faceI];


//     count_index = 0;
//     for (int i = 0; i < 3; i++)
//     {
//       for (int j = 0; j < 3; j++)
//       {

//         dsigmadt.boundaryFieldRef()[patchI][faceI].component(count_index) += sigAv(i, j); // sig(i,j);
//         D_P.boundaryFieldRef()[patchI][faceI].component(count_index) = D_P_Average(i, j); // F_D_P(i,j);
//         count_index++;
//       }
//     }



        
        
//         }
//     }
// }

  




  // for (label i = 0; i < N_Ori; i++)
  // {
  //   for (label countS = 0; countS < MatList[Phase2Material[i]].N_Slips; countS++)
  //   {
  //     g_alpha_dots[i][countS].correctBoundaryConditions();
  //     gamma_alpha_dots[i][countS].correctBoundaryConditions();
  //     // gamma_alphas[i][countS].correctBoundaryConditions();
  //     g_alphas[i][countS].correctBoundaryConditions();
  //   }
  // }



dsigmadt.correctBoundaryConditions();




}
