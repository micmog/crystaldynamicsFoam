    MRF.correctBoundaryVelocity(U);

    // fvc::interpolate(s);

    // surfaceTensorField sf("sf", fvc::interpolate(s));
    // volVectorField YT = fvc::div(sf);

    divs = fvc::div(s);

    // Scaling factor for Rhie-Chow stabilisation term
    // This is porportional to the material stiffness and time-step
    const dimensionedScalar k("k", dimPressure, 3e10);
    const scalar alpha(mesh.solution().lookupOrDefault<scalar>("alpha", 1.0));
    const dimensionedScalar stabScaleFactor(alpha*k*runTime.deltaT());

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      // + MRF.DDt(rho, U)
      //+ turbulence.divDevTau(rho, U)
      // - fvc::div(s)
         // simple stress model, compact stencil
      // - k*runTime.deltaT()*(fvc::laplacian(U) + fvc::div(gradU.T()))
         // less simple stress model, large stencil
        - k*runTime.deltaT()*(fvc::div(gradU + gradU.T() + tr(gradU)*I))
    ==
        stabScaleFactor*(fvm::laplacian(U) - fvc::div(gradU))
    //   phaseChange.SU(rho, rhoPhi, U)
    // + fvModels.source(rho, U)
    );

    UEqn.relax();

    fvConstraints.constrain(UEqn);

    UEqn.solve();

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    mixture.surfaceTensionForce()
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        fvConstraints.constrain(U);
    }

    U.relax();
