// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(U) + fvm::div(phi, U)
    + MRF.DDt(U)
    + turbulence->divDevReff(U)

 ==

    fvOptions(U)

);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

volScalarField UA(UEqn.A());
HFDIBDEM.interpolateIB(U,Ui,lambda);
Ui.correctBoundaryConditions();

f = surface*(UEqn.A()*Ui - UEqn.H() + fvc::grad(p));
f.correctBoundaryConditions();  

//~ if (pimple.momentumPredictor())
//~ {
    //~ solve(UEqn == -fvc::grad(p) + f);

    //~ fvOptions.correct(U);
//~ }

for (int i = 0; i < fDirectForcingIter; i++)
{
    Info << "Looping over UEqn, iteration: " << i << endl;
    solve(UEqn == -fvc::grad(p) + f);

    fvOptions.correct(U);
    
    HFDIBDEM.interpolateIB(U,Ui,lambda);
    
    if (max(mag(surface*(Ui - U))).value() < tolUEqn)
    {
        Info << "HFDIBDEM: U converged to Ui within max tolerance " << tolUEqn << endl;
        break;
    }
    
    Ui.correctBoundaryConditions();
    //~ f = surface*(UEqn.A()*Ui - UEqn.H() + fvc::grad(p));
    f += 0.8*surface*(Ui - U)/runTime.deltaT();                         //MI: in this loop, correct the force incrementaly
    f.correctBoundaryConditions(); 
    
    if (i == fDirectForcingIter - 1)
    {
        Info << "Loop over UEqn did not converge within " << fDirectForcingIter << " iterations" << endl;
    }
}

UEqn -= lambda*f;
//~ UEqn -= f;
