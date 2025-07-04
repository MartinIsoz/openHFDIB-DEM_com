volScalarField rAU(1.0/UEqn.A());
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::flux(HbyA)
  + MRF.zeroFilter(fvc::interpolate(rAU)*fvc::ddtCorr(U, phi, Uf))
);

MRF.makeRelative(phiHbyA);

if (p.needReference())
{
    fvc::makeRelative(phiHbyA, U);
    adjustPhi(phiHbyA, U, p);
    fvc::makeAbsolute(phiHbyA, U);
}

tmp<volScalarField> rAtU(rAU);

if (pimple.consistent())
{
    rAtU = 1.0/max(1.0/rAU - UEqn.H1(), 0.1/rAU);
    phiHbyA +=
        fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh.magSf();
    HbyA -= (rAU - rAtU())*fvc::grad(p);
}

if (pimple.nCorrPISO() <= 1)
{
    tUEqn.clear();
}

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, U, phiHbyA, rAtU(), MRF);


// Non-orthogonal pressure corrector loop
while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA)
    );

    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve();

    if (pimple.finalNonOrthogonalIter())
    {
        phi = phiHbyA - pEqn.flux();
    }
}

#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

//~ U = HbyA - rAtU*fvc::grad(p) + rAU*f;
U = HbyA - rAtU*fvc::grad(p);
U.correctBoundaryConditions();
fvOptions.correct(U);

// Correct Uf if the mesh is moving
fvc::correctUf(Uf, U, phi);

// Make the fluxes relative to the mesh motion
fvc::makeRelative(phi, U);

//~ fDragPress = fvc::grad(p);
//~ fDragVisc  = -fvc::div(turbulence->devReff());
//~ for (label pass=0; pass<=fDragSmoothingIter; pass++)
//~ {
    //~ fDragPress = fvc::average(fvc::interpolate(fDragPress));
    //~ fDragVisc  = fvc::average(fvc::interpolate(fDragVisc));
    //~ fDragPress.correctBoundaryConditions();
    //~ fDragVisc.correctBoundaryConditions();
//~ }

//~ HFDIBDEM.updateFSCoupling(lambda,fDragPress,fDragVisc);                 //MI: dubious line
HFDIBDEM.interpolateIB(U,Ui,lambda);

//Ui -= (Ui - U)*(1-lambda);

forAll(surface, sI)
{
    if (lambda[sI] > thrSurf)
        surface[sI] = 1;
    else
        surface[sI] = 0;
}

//~ f += 0.8*surface*(Ui - U)/runTime.deltaT();
//~ f = surface*(UA*Ui - UEqn.H() + fvc::grad(p));
//~ f.correctBoundaryConditions();
