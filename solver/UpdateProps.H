const volTensorField gradU(fvc::grad(U));
{

rhoPhi = fvc::interpolate(rho)*fvc::flux(U);
phi=fvc::flux(U);

dimensionedScalar shear_mod ("shear_mod", dimensionSet(1,-1,-2,0,0,0,0),80e9);

//volTensorField TgradU(gradU.T());
//TgradU.correctBoundaryConditions(); // already called by line above

// volTensorField TgradD = (fvc::grad(D)().T());
// TgradD.correctBoundaryConditions();
// gradU = fvc::grad(fvc::interpolate(U));


omega_dot = skew(gradU);
//omega_dot = 0.5*(fvc::grad(U)-TgradU); // this is the "skew" function, i.e. skew(gradU)
//omega_dot.correctBoundaryConditions(); // already called by line above
epsilon_dot = dev(symm(gradU));
// epsilon_dot = 0.5*(fvc::grad(U)+TgradU) - (1.0/3.0)*(fvc::div(phi))*I; // this looks like dev(gradU)
// epsilon_dot.correctBoundaryConditions(); // already called by line above

// epsilon = 0.5*(fvc::grad(D)-TgradD);
// epsilon.correctBoundaryConditions();
// dsdt = 2.0*shear_mod*(epsilon_dot/* - epsilon_P_dot*/) - ((s & omega_dot) - (omega_dot & s));
dsdt *=0.0;
dsdt = - ((s & omega_dot) - (omega_dot & s));


// surfaceTensorField epsilon_dotf = fvc::interpolate(epsilon_dot);

// C seems to be constant within this forAll loop; if so, this should be
// created once outside this loop
Fastor::Tensor<double,3,3,3,3> C(0.0);

double C11 = 8.0e10;
double C12 = 8.0e10;
double C44 = 8.0e10;


C(0,0,0,0)=C11;
C(1,1,1,1)=C11;
C(2,2,2,2)=C11;

C(1,0,1,0)=C44;
C(2,0,2,0)=C44;
C(0,0,1,1)=C12;
C(1,0,0,1)=C44;
C(0,0,2,2)=C12;
C(2,0,0,2)=C44;
C(0,1,1,0)=C44;
C(1,1,0,0)=C12;
C(0,1,0,1)=C44;
C(2,1,2,1)=C44;
C(1,1,2,2)=C12;
C(2,1,1,2)=C44;
C(0,2,2,0)=C44;
C(2,2,0,0)=C12;
C(1,2,2,1)=C44;
C(2,2,1,1)=C12;
C(0,2,0,2)=C44;
C(1,2,1,2)=C44;

dsdt_crystal*=0.0;
forAll(mesh.C(), cell)
{

// const labelList& faces = mesh.cells()[cell];

//     forAll( faces, facei )        // loop over all faces in cellID
//     {
//         // vector faceINormal = mesh.Sf()[facei] / mesh.magSf()[facei] ; 
//         // Info << " i = " << facei << ", faceINormal = " << faceINormal << endl ;
//         Info << " edot surface = " << facei << ": " << epsilon_dotf[facei] << endl ;





//     }

// Fastor::Tensor<double,3,3> I2={{1.0,0.0,0.0},
//                                  {0.0,1.0,0.0},
//                                  {0.0,0.0,1.0}};

// C = 8.0e10*Fastor::outer(I2,I2);


  // Fastor::print("C: ",C);
// exit(0);
//   Info<<"epsilon_dot: "<<epsilon_dot[cell]<<endl;

//  Info<<"epsilon_dot 4: "<<epsilon_dot[cell].component(4)<<endl;
//  Info<<"epsilon_dot 5: "<<epsilon_dot[cell].component(5)<<endl;

  Fastor::Tensor<double,3,3> Fepsilon_dot(0.0);
  int count_index=0;
  for(int i = 0 ; i < 3 ; i++ ){
    for(int j = 0 ; j < 3 ; j++ ){

      Fepsilon_dot(i,j)=epsilon_dot[cell].component(count_index);
      // Fepsilon_dot(i,j)=epsilon[cell].component(count_index);

      count_index++;
    }
  }



// Fastor::Tensor<double,9,9> reshapedC(0.0) ;
// reshapedC = Fastor::reshape<9,9>(C);

// Fastor::print("reshapedC: ",reshapedC);

// exit(0);

enum {I,J,K,L};
Fastor::Tensor<double,3,3> sig= Fastor::einsum<Fastor::Index<I,J,K,L>,Fastor::Index<K,L>/*,Fastor::OIndex<I,J>*/>(C,Fepsilon_dot);////Double contraction


// if(mag(epsilon_dot[cell])>SMALL){

// Fastor::print("edot ",Fepsilon_dot);
// Fastor::print("sig",sig);

// }


count_index=0;
for(int i = 0 ; i < 3 ; i++ ){
    for(int j = 0 ; j < 3 ; j++ ){

      dsdt_crystal[cell].component(count_index) =sig(i,j);
      // testSIG[cell].component(count_index) =sig(i,j);
  count_index++;
    }

}


}
dsdt_crystal.correctBoundaryConditions();
dsdt+=dsdt_crystal;
//dsdt = 1.0*dsdt_crystal; // no Jaumann rate correction

// dsdt_crystal.write();
// exit(0);


dsdt.correctBoundaryConditions();
s=s + runTime.deltaT()*(dsdt);

    // solve
    // (
    //     fvm::ddt(s) ==
    //     dsdt
    // );

s.correctBoundaryConditions();




  //   fvVectorMatrix DEqn
  //   (
  //       fvm::d2dt2(rho, D) 
  //     -fvc::div(testSIG)
  //  //  ==
  //   //   phaseChange.SU(rho, rhoPhi, U)
  //   // + fvModels.source(rho, U)
  //   );

  //   // DEqn.relax();

  //   // fvConstraints.constrain(DEqn);

  //   DEqn.solve();



// for (label count = 0; count < N_Ori; count++)
// {}


//TESTING AREA

// Fastor::Tensor<double,3,3> teststrain={{1.0,2.0,3.0},
//                                  {4.0,5.0,6.0},
//                                  {7.0,8.0,9.0}};

// Fastor::Tensor<double,3,3> sigdiv={{2.0,4.0,6.0},
//                                  {8.0,10.0,12.0},
//                                  {14.0,16.0,18.0}};





// Fastor::Tensor<double,3,3,3,3> Stiffness_matrix(0.0);

//   // Stiffness_matrix(0,0,0,0)=2.0;
//   // Stiffness_matrix(0,0,1,1)=0.0;
//   // Stiffness_matrix(0,0,2,2)=0.0;
//   // Stiffness_matrix(1,1,0,0)=0.0;
//   // Stiffness_matrix(1,1,1,1)=2.0;
//   // Stiffness_matrix(1,1,2,2)=0.0;
//   // Stiffness_matrix(2,2,0,0)=0.0;
//   // Stiffness_matrix(2,2,1,1)=0.0;
//   // Stiffness_matrix(2,2,2,2)=2.0;

//   Stiffness_matrix(0,0,0,0)=2.0;
//   Stiffness_matrix(0,0,0,1)=3.0;
//   Stiffness_matrix(0,0,0,2)=4.0;
//   Stiffness_matrix(0,0,1,0)=5.0;
//   Stiffness_matrix(0,0,1,1)=6.0;
//   Stiffness_matrix(0,0,1,2)=7.0;
//   Stiffness_matrix(0,0,2,0)=8.0;
//   Stiffness_matrix(0,0,2,1)=9.0;
//   Stiffness_matrix(0,0,2,2)=10.0;

//   Stiffness_matrix(0,1,0,0)=12.0;
//   Stiffness_matrix(0,1,0,1)=13.0;
//   Stiffness_matrix(0,1,0,2)=14.0;
//   Stiffness_matrix(0,1,1,0)=15.0;
//   Stiffness_matrix(0,1,1,1)=16.0;
//   Stiffness_matrix(0,1,1,2)=17.0;
//   Stiffness_matrix(0,1,2,0)=18.0;
//   Stiffness_matrix(0,1,2,1)=19.0;
//   Stiffness_matrix(0,1,2,2)=110.0;

//   Stiffness_matrix(2,2,0,0)=12.0;
//   Stiffness_matrix(2,2,0,1)=13.0;
//   Stiffness_matrix(2,2,0,2)=14.0;
//   Stiffness_matrix(2,2,1,0)=15.0;
//   Stiffness_matrix(2,2,1,1)=16.0;
//   Stiffness_matrix(2,2,1,2)=17.0;
//   Stiffness_matrix(2,2,2,0)=18.0;
//   Stiffness_matrix(2,2,2,1)=19.0;
//   Stiffness_matrix(2,2,2,2)=110.0;


// enum {I,J,K,L};
// Fastor::Tensor<double,3,3> TEST_OUT= Fastor::einsum<Fastor::Index<I,J,K,L>,Fastor::Index<K,L>/*,Fastor::OIndex<I,J>*/>(Stiffness_matrix,teststrain);//

// Fastor::print("TEST_stress", TEST_OUT);

// exit(0);





// for(int i = 0 ; i < 3 ; i++ ){
//   for(int j = 0 ; j < 3 ; j++ ){
//     for(int k = 0 ; k < 3 ; k++ ){
//       for(int l = 0 ;l  < 3 ; l++ ){

//         double value =0.0;
//         for(int m = 0 ; m < 3 ; m++ ){
//           for(int n = 0 ; n < 3 ; n++ ){
//             for(int p = 0 ; p < 3 ; p++ ){
//               for(int q = 0 ;q  < 3 ; q++ ){

//                     value+=Rot_matrix(i,m)*Rot_matrix(j,n)*Rot_matrix(k,p)*Rot_matrix(l,q)*C(m,n,p,q);

//                     }
//                   }
//                 }
//               }

//         Cprime(i,j,k,l)=value;
      
//       }
//     }
//   }
// }

  // exit(0);

}
