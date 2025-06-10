// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(U) + fvm::div(phi, U)
  + MRF.DDt(U)

//zlobi divDevSigma
//  + turbulence->divDevSigma(U)

//nezlobi
    + turbulence->divDevReff(U)

 ==

    fvOptions(U)

);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p) + f);

    fvOptions.correct(U);
}
