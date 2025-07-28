{


Euler_Angles*=0.0;

  forAll(mesh.C(), cell)
  {

        label largest_field=-1;
    double current_largest_value = -1.0;

    for (label h = 0; h < N_Ori; h++)
    {

      if(PopBal[h][cell]>current_largest_value){
        current_largest_value = PopBal[h][cell];
      
        largest_field = h;
      }


            double phi_1_total = Rot1[h][cell]; // ROT_1[cell];// phi1[h] ;//+ omega_tot[cell].component(5);

            double phi_2_total = Rot2[h][cell]; // ROT_2[cell];// phi2[h] ;//+ omega_tot[cell].component(2);

            double phi_3_total = Rot3[h][cell]; // ROT_3[cell];// phi3[h] ;//+ omega_tot[cell].component(1);


            Fastor::Tensor<double, 3, 3> Rot_matrix(0.0);

      // Rot_matrix(0, 0) = test_rot_diff[h][cell].component(0);// Foam::cos(phi_2_total) * Foam::cos(phi_3_total);
      // Rot_matrix(0, 1) = test_rot_diff[h][cell].component(1);//(Foam::sin(phi_1_total) * Foam::sin(phi_2_total) * Foam::cos(phi_3_total)) - (Foam::cos(phi_1_total) * Foam::sin(phi_3_total));
      // Rot_matrix(0, 2) = test_rot_diff[h][cell].component(2);//(Foam::cos(phi_1_total) * Foam::sin(phi_2_total) * Foam::cos(phi_3_total)) + (Foam::sin(phi_1_total) * Foam::sin(phi_3_total));

      // Rot_matrix(1, 0) = test_rot_diff[h][cell].component(3);//Foam::cos(phi_2_total) * Foam::sin(phi_3_total);
      // Rot_matrix(1, 1) = test_rot_diff[h][cell].component(4);//(Foam::sin(phi_1_total) * Foam::sin(phi_2_total) * Foam::sin(phi_3_total)) + (Foam::cos(phi_1_total) * Foam::cos(phi_3_total));
      // Rot_matrix(1, 2) = test_rot_diff[h][cell].component(5);//(Foam::cos(phi_1_total) * Foam::sin(phi_2_total) * Foam::sin(phi_3_total)) - (Foam::sin(phi_1_total) * Foam::cos(phi_3_total));

      // Rot_matrix(2, 0) = test_rot_diff[h][cell].component(6);//-Foam::sin(phi_2_total);
      // Rot_matrix(2, 1) =test_rot_diff[h][cell].component(7);// Foam::sin(phi_1_total) * Foam::cos(phi_2_total);
      // Rot_matrix(2, 2) = test_rot_diff[h][cell].component(8);//Foam::cos(phi_1_total) * Foam::cos(phi_2_total);

      Rot_matrix(0, 0) = Foam::cos(phi_2_total) * Foam::cos(phi_3_total);
      Rot_matrix(0, 1) = (Foam::sin(phi_1_total) * Foam::sin(phi_2_total) * Foam::cos(phi_3_total)) - (Foam::cos(phi_1_total) * Foam::sin(phi_3_total));
      Rot_matrix(0, 2) = (Foam::cos(phi_1_total) * Foam::sin(phi_2_total) * Foam::cos(phi_3_total)) + (Foam::sin(phi_1_total) * Foam::sin(phi_3_total));

      Rot_matrix(1, 0) = Foam::cos(phi_2_total) * Foam::sin(phi_3_total);
      Rot_matrix(1, 1) = (Foam::sin(phi_1_total) * Foam::sin(phi_2_total) * Foam::sin(phi_3_total)) + (Foam::cos(phi_1_total) * Foam::cos(phi_3_total));
      Rot_matrix(1, 2) = (Foam::cos(phi_1_total) * Foam::sin(phi_2_total) * Foam::sin(phi_3_total)) - (Foam::sin(phi_1_total) * Foam::cos(phi_3_total));

      Rot_matrix(2, 0) = -Foam::sin(phi_2_total);
      Rot_matrix(2, 1) = Foam::sin(phi_1_total) * Foam::cos(phi_2_total);
      Rot_matrix(2, 2) = Foam::cos(phi_1_total) * Foam::cos(phi_2_total);



      // Fastor::print("Rot Matrix: ",Rot_matrix);

      // https://en.wikipedia.org/wiki/Rotation_matrix -
      //"represents an extrinsic rotation whose (improper) Euler angles are α, β, γ, about axes x, y, z" - phi1 ==alpha, phi2==beta,phi3==gamma

        double deltaT_ = /*(1.0 / n_sub_steps) */ runTime.deltaT().value(); 
        // Info<<deltaT_<<endl;
        Fastor::Tensor<double, 3, 3> Elastic_Spin_h_deltaT(0.0);
        Fastor::Tensor<double, 3, 3> Elastic_Spin_h_deltaT_sq(0.0);
        int count_index_E = 0;
        for (int ii = 0; ii < 3; ii++)
        {
          for (int jj = 0; jj < 3; jj++)
          {
        // Info<<i<<","<<j<<endl;
            Elastic_Spin_h_deltaT(ii, jj) = (omega_dot[cell].component(count_index_E) - Omega_P[h][cell].component(count_index_E));//shouldnt be the elastic spin, should be the total spin elastic + plastic

            count_index_E++;
          }
        }
        // Info<<omega_dot[cell].component(7)<<endl;;
        // Fastor::print("TEST",Elastic_Spin_h_deltaT);
        Elastic_Spin_h_deltaT_sq = Fastor::matmul(Elastic_Spin_h_deltaT,Elastic_Spin_h_deltaT);  



     //Update rotation Matrix Based on Elastic spin

     Fastor::Tensor<double, 3, 3> Rot_matrix_N_plus_one(0.0);
     Fastor::Tensor<double, 3, 3> I_d(0.0);
     for(int set=0;set<3;set++){
      I_d(set,set)=1.0;
     }
      
    //  Fastor::print("I",I_d);

    double X_comps = Foam::sqrt(
                    // Foam::pow(Elastic_Spin_h_deltaT(0,0),2.0)
                    // +Foam::pow(Elastic_Spin_h_deltaT(0, 1),2.0)
                    // +Foam::pow(Elastic_Spin_h_deltaT(0, 2),2.0)
                    // +Foam::pow(Elastic_Spin_h_deltaT(1, 0),2.0)
                    // +Foam::pow(Elastic_Spin_h_deltaT(1, 1),2.0)
                    // +Foam::pow(Elastic_Spin_h_deltaT(1, 2),2.0)
                    // +Foam::pow(Elastic_Spin_h_deltaT(2, 0),2.0)
                    // +Foam::pow(Elastic_Spin_h_deltaT(2, 1),2.0)
                    // +Foam::pow(Elastic_Spin_h_deltaT(2,2),2.0)

                    Foam::pow(Elastic_Spin_h_deltaT(2,1),2.0)
                    +Foam::pow(Elastic_Spin_h_deltaT(0, 2),2.0)
                    +Foam::pow(Elastic_Spin_h_deltaT(1,0),2.0)
                    
                    );

    //                 Fastor::print("Elastic_Spin_h_deltaT",Elastic_Spin_h_deltaT);

    // Info<<"X_comps: "<<X_comps<<endl;
    Fastor::Tensor<double, 3, 3> e_to_the_C(0.0);
    e_to_the_C = I_d;
    if(X_comps>SMALL){
      // Info<<"x_comps"<<X_comps<<endl;
    

    e_to_the_C  +=  ((Foam::sin(X_comps*deltaT_)/X_comps)*Elastic_Spin_h_deltaT) + ((1.0-Foam::cos(X_comps*deltaT_))/Foam::pow(X_comps,2.0))* Elastic_Spin_h_deltaT_sq;


    
    }

    //this below  is for the Openfoam implementation
    // Rot_matrix_N_plus_one = Fastor::matmul(Rot_matrix,Fastor::transpose(e_to_the_C));//If we define euler angles relative to +ve rotation around axes then this update has to be consistent with this
    //e_to_the_C is the rotation update, now need to rotate the rotation matrix by this

  Rot_matrix_N_plus_one = Fastor::matmul(e_to_the_C,Rot_matrix);

    // test_rot_diff[h][cell].component(0) = Rot_matrix_N_plus_one(0,0);
    // test_rot_diff[h][cell].component(1) = Rot_matrix_N_plus_one(0,1);
    // test_rot_diff[h][cell].component(2) = Rot_matrix_N_plus_one(0,2);
    // test_rot_diff[h][cell].component(3) = Rot_matrix_N_plus_one(1,0);
    // test_rot_diff[h][cell].component(4) = Rot_matrix_N_plus_one(1,1);
    // test_rot_diff[h][cell].component(5) = Rot_matrix_N_plus_one(1,2);
    // test_rot_diff[h][cell].component(6) = Rot_matrix_N_plus_one(2,0);
    // test_rot_diff[h][cell].component(7) = Rot_matrix_N_plus_one(2,1);
    // test_rot_diff[h][cell].component(8) = Rot_matrix_N_plus_one(2,2);
    
    // if(cell==240&&h==1){
    //   Info<<h<<endl;
    // Fastor::print("Rot: ",Rot_matrix);
    // Fastor::print("e_to_the_C: ",e_to_the_C);

    // Fastor::print("Rot_n+1: ",Rot_matrix_N_plus_one);
    // // exit(0);
    // }

    // Rot1[h][cell] = ?
    // Rot2[h][cell] = ?
    // Rot3[h][cell] = ?

    double theta_angle = 0.0;
    double psi_angle = 0.0;
    double phi_angle = 0.0;

    if(mag(Rot_matrix_N_plus_one(2,0))!=1.0){
      // Info<<"HERE"<<endl;
      double theta_1 = -std::asin(Rot_matrix_N_plus_one(2,0));
      double theta_2 = M_PI - theta_1;
      double psi_1 = std::atan2((Rot_matrix_N_plus_one(2,1)/std::cos(theta_1)),(Rot_matrix_N_plus_one(2,2)/std::cos(theta_1)));
      double psi_2 = std::atan2((Rot_matrix_N_plus_one(2,1)/std::cos(theta_2)),(Rot_matrix_N_plus_one(2,2)/std::cos(theta_2)));
      double phi_1 = std::atan2((Rot_matrix_N_plus_one(1,0)/std::cos(theta_1)),(Rot_matrix_N_plus_one(0,0)/std::cos(theta_1)));
      double phi_2 = std::atan2((Rot_matrix_N_plus_one(1,0)/std::cos(theta_2)),(Rot_matrix_N_plus_one(0,0)/std::cos(theta_2)));

      // if(cell==240/*&&h==1*/){
      // Info<<"angles "<<h<<" "<<endl;
      // Info<<"theta1: "<<theta_1<<endl;
      // Info<<"theta_2: "<<theta_2<<endl;
      // Info<<"psi_1: "<<psi_1<<endl;
      // Info<<"psi_2: "<<psi_2<<endl;
      // Info<<"phi_1: "<<phi_1<<endl;
      // Info<<"phi_2: "<<phi_2<<endl;
      // Info<<""<<endl;
      // }

      theta_angle = theta_1;//max(theta_1,theta_2);
      psi_angle = psi_1;//max(psi_1,psi_2);
      phi_angle = phi_1;//max(phi_1,phi_2);

      // Info<<"TEST: "<<psi_1<<"\t"<<psi_2<<endl;

    }
    else{
      phi_angle = 0.0;
        if(Rot_matrix_N_plus_one(2,0)==-1.0){
          theta_angle = M_PI/2.0;
          psi_angle = phi_angle + std::atan2(Rot_matrix_N_plus_one(0,1),(Rot_matrix_N_plus_one(0,2)));
        }
        else{
          theta_angle = -M_PI/2.0;
          psi_angle = -phi_angle + std::atan2(-Rot_matrix_N_plus_one(0,1),-Rot_matrix_N_plus_one(0,2));
        }

    }

    Rot1[h][cell] = psi_angle;//update euler angles for advection
    Rot2[h][cell] = theta_angle;
    Rot3[h][cell] = phi_angle;




    }

          Euler_Angles[cell].component(0) = Rot1[largest_field][cell];
      Euler_Angles[cell].component(1) = Rot2[largest_field][cell];
      Euler_Angles[cell].component(2) = Rot3[largest_field][cell];

  }








//Over Boundaries
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


                    label largest_field=-1;
    double current_largest_value = -1.0;


    for (label h = 0; h < N_Ori; h++)
    {


        if(PopBal[h].boundaryField()[patchI][faceI]>current_largest_value){
        current_largest_value = PopBal[h].boundaryField()[patchI][faceI];
      
        largest_field = h;
      }


          double phi_1_total = Rot1[h].boundaryField()[patchI][faceI]; // ROT_1[cell];// phi1[h] ;//+ omega_tot[cell].component(5);

      double phi_2_total = Rot2[h].boundaryField()[patchI][faceI]; // ROT_2[cell];// phi2[h] ;//+ omega_tot[cell].component(2);

      double phi_3_total = Rot3[h].boundaryField()[patchI][faceI]; // ROT_3[cell];// phi3[h] ;//+ omega_tot[cell].component(1);


      Fastor::Tensor<double, 3, 3> Rot_matrix(0.0);

      Rot_matrix(0, 0) = Foam::cos(phi_2_total) * Foam::cos(phi_3_total);
      Rot_matrix(0, 1) = (Foam::sin(phi_1_total) * Foam::sin(phi_2_total) * Foam::cos(phi_3_total)) - (Foam::cos(phi_1_total) * Foam::sin(phi_3_total));
      Rot_matrix(0, 2) = (Foam::cos(phi_1_total) * Foam::sin(phi_2_total) * Foam::cos(phi_3_total)) + (Foam::sin(phi_1_total) * Foam::sin(phi_3_total));

      Rot_matrix(1, 0) = Foam::cos(phi_2_total) * Foam::sin(phi_3_total);
      Rot_matrix(1, 1) = (Foam::sin(phi_1_total) * Foam::sin(phi_2_total) * Foam::sin(phi_3_total)) + (Foam::cos(phi_1_total) * Foam::cos(phi_3_total));
      Rot_matrix(1, 2) = (Foam::cos(phi_1_total) * Foam::sin(phi_2_total) * Foam::sin(phi_3_total)) - (Foam::sin(phi_1_total) * Foam::cos(phi_3_total));

      Rot_matrix(2, 0) = -Foam::sin(phi_2_total);
      Rot_matrix(2, 1) = Foam::sin(phi_1_total) * Foam::cos(phi_2_total);
      Rot_matrix(2, 2) = Foam::cos(phi_1_total) * Foam::cos(phi_2_total);
      // https://en.wikipedia.org/wiki/Rotation_matrix -
      //"represents an extrinsic rotation whose (improper) Euler angles are α, β, γ, about axes x, y, z" - phi1 ==alpha, phi2==beta,phi3==gamma


        double deltaT_ = /*(1.0 / n_sub_steps) */ runTime.deltaT().value(); 

        Fastor::Tensor<double, 3, 3> Elastic_Spin_h_deltaT(0.0);
        Fastor::Tensor<double, 3, 3> Elastic_Spin_h_deltaT_sq(0.0);
        int count_index_E = 0;
        for (int ii = 0; ii < 3; ii++)
        {
          for (int jj = 0; jj < 3; jj++)
          {
        // Info<<i<<","<<j<<endl;
            Elastic_Spin_h_deltaT(ii, jj) = deltaT_*(omega_dot.boundaryField()[patchI][faceI].component(count_index_E) - Omega_P[h].boundaryField()[patchI][faceI].component(count_index_E));//
            // Elastic_Spin_h_sq(ii, jj) = Foam::pow(omega_dot[cell].component(count_index_E) - Omega_P[h][cell].component(count_index_E),2.0);
            count_index_E++;
          }
        }

        Elastic_Spin_h_deltaT_sq = Fastor::matmul(Elastic_Spin_h_deltaT,Elastic_Spin_h_deltaT);  



     //Update rotation Matrix Based on Elastic spin

     Fastor::Tensor<double, 3, 3> Rot_matrix_N_plus_one(0.0);
     Fastor::Tensor<double, 3, 3> I_d(0.0);
     for(int set=0;set<3;set++){
      I_d(set,set)=1.0;
     }
      
    //  Fastor::print("I",I_d);

    double X_comps = Foam::sqrt(Foam::pow(Elastic_Spin_h_deltaT(2,1),2.0)
                    +Foam::pow(Elastic_Spin_h_deltaT(0, 2),2.0)
                    +Foam::pow(Elastic_Spin_h_deltaT(1,0),2.0));

                    // Fastor::print("Elastic_Spin_h_deltaT",Elastic_Spin_h_deltaT);

    // Info<<"X_comps: "<<X_comps<<endl;
    Fastor::Tensor<double, 3, 3> e_to_the_C(0.0);
    e_to_the_C = I_d;
    if(X_comps>SMALL){
      // Info<<"x_comps"<<X_comps<<endl;
    

    e_to_the_C  += ((Foam::sin(X_comps)/X_comps)*Elastic_Spin_h_deltaT) + ((1.0-Foam::cos(X_comps))/Foam::pow(X_comps,2.0))* Elastic_Spin_h_deltaT_sq;


    
    }

    //this below  is for the Openfoam implementation
    // Rot_matrix_N_plus_one = Fastor::matmul(Rot_matrix,Fastor::transpose(e_to_the_C));//If we define euler angles relative to +ve rotation around axes then this update has to be consistent with this
    Rot_matrix_N_plus_one = Fastor::matmul(e_to_the_C,Rot_matrix);


    // Rot1[h][cell] = ?
    // Rot2[h][cell] = ?
    // Rot3[h][cell] = ?

    double theta_angle = 0.0;
    double psi_angle = 0.0;
    double phi_angle = 0.0;

    if(mag(Rot_matrix_N_plus_one(2,0))!=1.0){
      double theta_1 = -std::asin(Rot_matrix_N_plus_one(2,0));
      double theta_2 = M_PI - theta_1;
      double psi_1 = std::atan2((Rot_matrix_N_plus_one(2,1)/std::cos(theta_1)),(Rot_matrix_N_plus_one(2,2)/std::cos(theta_1)));
      double psi_2 = std::atan2((Rot_matrix_N_plus_one(2,1)/std::cos(theta_2)),(Rot_matrix_N_plus_one(2,2)/std::cos(theta_2)));
      double phi_1 = std::atan2((Rot_matrix_N_plus_one(1,0)/std::cos(theta_1)),(Rot_matrix_N_plus_one(0,0)/std::cos(theta_1)));
      double phi_2 = std::atan2((Rot_matrix_N_plus_one(1,0)/std::cos(theta_2)),(Rot_matrix_N_plus_one(0,0)/std::cos(theta_2)));

      theta_angle = theta_1;//max(theta_1,theta_2);
      psi_angle = psi_1;//max(psi_1,psi_2);
      phi_angle = phi_1;//max(phi_1,phi_2);

      // Info<<"TEST: "<<psi_1<<"\t"<<psi_2<<endl;

    }
    else{
      phi_angle = 0.0;
        if(Rot_matrix_N_plus_one(2,0)==-1.0){
          theta_angle = M_PI/2.0;
          psi_angle = phi_angle + std::atan2(Rot_matrix_N_plus_one(0,1),(Rot_matrix_N_plus_one(0,2)));
        }
        else{
          theta_angle = -M_PI/2.0;
          psi_angle = -phi_angle + std::atan2(-Rot_matrix_N_plus_one(0,1),-Rot_matrix_N_plus_one(0,2));
        }

    }

    Rot1[h].boundaryFieldRef()[patchI][faceI] = psi_angle;//update euler angles for advection
    Rot2[h].boundaryFieldRef()[patchI][faceI] = theta_angle;
    Rot3[h].boundaryFieldRef()[patchI][faceI] = phi_angle;





    }//END LOOP OVER FIELDS

          Euler_Angles.boundaryFieldRef()[patchI][faceI].component(0) = Rot1[largest_field].boundaryField()[patchI][faceI];
      Euler_Angles.boundaryFieldRef()[patchI][faceI].component(1) = Rot2[largest_field].boundaryField()[patchI][faceI];
      Euler_Angles.boundaryFieldRef()[patchI][faceI].component(2) = Rot3[largest_field].boundaryField()[patchI][faceI];








        }
    }
}//END OVER BOUNDARIES

  



}