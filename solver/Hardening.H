
    for (label i = 0; i < N_Ori; i++)
  {
    for (label countS = 0; countS < MatList[Phase2Material[i]].N_Slips; countS++)
    {
      g_alpha_dots[i][countS] *= 0.0;
      ;
      // gamma_alpha_dots[i][countS].correctBoundaryConditions();
      // gamma_alphas[i][countS].correctBoundaryConditions();
      // g_alphas[i][countS].correctBoundaryConditions();
    }
  }



// // // //   // // //update hardening rates
// // // // // // // Info<<"HERE1"<<endl;

//   forAll(mesh.C(), cell)
//   {
//     //   label Mat_name =-1;
//     //   double max_phase = -1.0;
//     for (label h = 0; h < N_Ori; h++)
//     {
//       if (PopBal[h][cell] > SMALL)
//       {

//         for (int S_alpha = 0; S_alpha < MatList[Phase2Material[h]].N_Slips; S_alpha++)
//         {


//           g_alpha_dots[h][S_alpha][cell] *= 0.0;

//           for (int S_beta = 0; S_beta < MatList[Phase2Material[h]].N_Slips; S_beta++)
//           {

//             //   // double h_of_gamma = ho*Foam::pow(1.0/(Foam::cosh((ho*Acc_Strain[cell])/(ta_s-tau_0))),2.0);
//             double h_alpha_beta = 0.0;
//             // Pout<<Pstream::myProcNo()<<"test_strain: "<<Acc_Strain[cell]<<endl;

//             double argument = (MatList[Phase2Material[h]].h_0 * Acc_Strain[cell]/*gamma_alphas[h][S_alpha][cell]*/) / (MatList[Phase2Material[h]].tau_S - MatList[Phase2Material[h]].initial_CRSS);
//             // argument = max(min(argument,20.0),-20.0);

//             // double h_of_gamma = MatList[Phase2Material[h]].h_0 * Foam::pow(1.0 / (Foam::cosh((MatList[Phase2Material[h]].h_0 * gamma_alphas[h][S_alpha][cell]) / (MatList[Phase2Material[h]].tau_S - MatList[Phase2Material[h]].initial_CRSS))), 2.0);
//               double h_of_gamma = MatList[Phase2Material[h]].h_0 * Foam::pow(1.0 / (Foam::cosh(argument)), 2.0);            

//             if (S_beta == S_alpha)
//             {
//               h_alpha_beta = h_of_gamma * (MatList[Phase2Material[h]].r + (1.0 - MatList[Phase2Material[h]].r));
//             }
//             else
//             {
//               h_alpha_beta = h_of_gamma * (MatList[Phase2Material[h]].r);
//             }

//             g_alpha_dots[h][S_alpha][cell] += /*PopBal[h][cell]*/(h_alpha_beta)*gamma_alpha_dots[h][S_beta][cell];
//           }
//           // double deltaT_ = /*(1.0 / n_sub_steps) */ runTime.deltaT().value(); 
//           // // g_alpha_dots[h][S_alpha][cell]=g_alpha_dots[h][S_alpha][cell]/deltaT_;
//           // g_alpha_dots[h][S_alpha][cell]=max(min(g_alpha_dots[h][S_alpha][cell],0.99/deltaT_),-0.99/deltaT_);
//           //
//         }

//     }

//     }
//   }


// //Over Boundaries
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


//    for (label h = 0; h < N_Ori; h++)
//     {
//             if (PopBal[h].boundaryField()[patchI][faceI] > SMALL)
//       {


//         for (int S_alpha = 0; S_alpha < MatList[Phase2Material[h]].N_Slips; S_alpha++)
//         {


//           g_alpha_dots[h][S_alpha].boundaryFieldRef()[patchI][faceI] *= 0.0;

//           for (int S_beta = 0; S_beta < MatList[Phase2Material[h]].N_Slips; S_beta++)
//           {


//             double h_alpha_beta = 0.0;

//             double argument = (MatList[Phase2Material[h]].h_0 * Acc_Strain.boundaryField()[patchI][faceI]/*gamma_alphas[h][S_alpha].boundaryField()[patchI][faceI]*/) / (MatList[Phase2Material[h]].tau_S - MatList[Phase2Material[h]].initial_CRSS);
//             // argument = max(min(argument,20.0),-20.0);

//               double h_of_gamma = MatList[Phase2Material[h]].h_0 * Foam::pow(1.0 / (Foam::cosh(argument)), 2.0);            

//             if (S_beta == S_alpha)
//             {
//               h_alpha_beta = h_of_gamma * (MatList[Phase2Material[h]].r + (1.0 - MatList[Phase2Material[h]].r));
//             }
//             else
//             {
//               h_alpha_beta = h_of_gamma * (MatList[Phase2Material[h]].r);
//             }

//             g_alpha_dots[h][S_alpha].boundaryFieldRef()[patchI][faceI] += /*PopBal[h].boundaryField()[patchI][faceI]*/(h_alpha_beta)*gamma_alpha_dots[h][S_beta].boundaryField()[patchI][faceI];
//           }
//           // double deltaT_ = /*(1.0 / n_sub_steps) */ runTime.deltaT().value(); 
//           // // g_alpha_dots[h][S_alpha].boundaryFieldRef()[patchI][faceI]=g_alpha_dots[h][S_alpha].boundaryFieldRef()[patchI][faceI]/deltaT_;
//           // g_alpha_dots[h][S_alpha].boundaryFieldRef()[patchI][faceI]=max(min(g_alpha_dots[h][S_alpha].boundaryField()[patchI][faceI],0.99/deltaT_),-0.99/deltaT_);
//           //
//         }

//     }


//     }



//         }

//     }

// }


// // //END update hardening rates


