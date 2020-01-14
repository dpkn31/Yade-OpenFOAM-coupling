/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    DPMFoam

Description
    Transient solver for the coupled transport of a single kinematic particle
    cloud including the effect of the volume fraction of particles on the
    continuous phase.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "PhaseIncompressibleTurbulenceModel.H"
#include "pimpleControl.H"
#include "FoamYade.H"

int main(int argc, char *argv[])
{
  
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    bool gaussianInterp = true;
    FoamYade yadeCoupling(mesh,Uc, gradP, vGrad, divT,ddtU_f,g,uSourceDrag,alphac,uSource,uParticle,gaussianInterp);
    yadeCoupling.setScalarProperties(partDensity.value(), rhocValue.value(), nuValue.value());
    
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;
	
        continuousPhaseTransport.correct();
        muc = rhoc*continuousPhaseTransport.nu();

	ddtU_f = fvc::ddt(Uc)+fvc::div(phic, Uc);
        gradP = fvc::grad(p);
        divT = 2*nuValue.value()*fvc::laplacian(alphac, Uc);
        vGrad = fvc::grad(Uc);

        yadeCoupling.setParticleAction(runTime.deltaT().value());



        // Update continuous phase volume fraction field
        alphac.correctBoundaryConditions();
        alphacf = fvc::interpolate(alphac);
	alphaPhic = alphacf*phic;
	
	uSource.correctBoundaryConditions(); 
	uSourceDrag.correctBoundaryConditions(); 

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UcEqn.H"

            // --- PISO loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                continuousPhaseTurbulence->correct();
            }
        }

        runTime.write();
	
	yadeCoupling.setSourceZero(); 
	
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
