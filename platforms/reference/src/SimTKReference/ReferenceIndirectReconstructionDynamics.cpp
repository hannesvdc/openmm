/*
 * ReferenceIndirectReconstructionDynamics.cpp
 *
 *  Created on: 2 Nov 2020
 *      Author: hannesvdc
 */

#include <cstring>
#include <sstream>
#include <iostream>

#include "SimTKOpenMMUtilities.h"
#include "ReferenceIndirectReconstructionDynamics.h"
#include "ReferenceVirtualSites.h"
#include "openmm/OpenMMException.h"

#include <cstdio>

using std::vector;
using namespace OpenMM;

/**---------------------------------------------------------------------------------------

   ReferenceBrownianDynamics constructor

   @param numberOfAtoms  number of atoms
   @param deltaT         delta t for dynamics
   @param temperature    temperature

   --------------------------------------------------------------------------------------- */

ReferenceIndirectReconstructionDynamics::ReferenceIndirectReconstructionDynamics(int numberOfAtoms,
                                                                                 double deltaT,
                                                                                 double temperature,
																			   ReactionCoordinate* rc,
																			   double lam ) :
           ReferenceDynamics(numberOfAtoms, deltaT, temperature), lambda(lam), reactionCoordinate(rc) {

   xPrime.resize(numberOfAtoms);
    macroVariable.resize(numberOfAtoms, Vec3(0., 0., 0.));
}

/**---------------------------------------------------------------------------------------

   ReferenceBrownianDynamics destructor

   --------------------------------------------------------------------------------------- */

ReferenceIndirectReconstructionDynamics::~ReferenceIndirectReconstructionDynamics() {
}


/**---------------------------------------------------------------------------------------

   Update -- driver routine for performing Brownian dynamics update of coordinates
   and velocities

   @param system              the System to be integrated
   @param atomCoordinates     atom coordinates
   @param velocities          velocities
   @param forces              forces
   @param masses              atom masses

   --------------------------------------------------------------------------------------- */

void ReferenceIndirectReconstructionDynamics::update(const OpenMM::System& system, vector<Vec3>& atomCoordinates,
                                                    vector<Vec3>& velocities,
                                                    vector<Vec3>& forces, vector<double>& masses, double tolerance) {

   // Perform the integration.
   int numberOfAtoms = system.getNumParticles();
    double beta = sqrt(2.0*BOLTZ*getTemperature());
   const double noiseAmplitude = beta*sqrt(getDeltaT());

   std::vector<Vec3> rcValue = reactionCoordinate->value(atomCoordinates);
   for (int i = 0; i < macroVariable.size(); ++i) {
	   rcValue[i][0] -= macroVariable[i][0];
	   rcValue[i][1] -= macroVariable[i][1];
	   rcValue[i][2] -= macroVariable[i][2];
   }
   std::vector<Vec3> rcGrad = reactionCoordinate->gradMatMul(atomCoordinates, rcValue);

   for (int i = 0; i < numberOfAtoms; ++i) {
       if (masses[i] != 0.0)
           for (int j = 0; j < 3; ++j) {
               xPrime[i][j] = atomCoordinates[i][j] + getDeltaT()*forces[i][j] - getDeltaT()*getLambda()*rcGrad[i][j] + noiseAmplitude*SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();
           }
   }

   // Update the positions and velocities.
   double velocityScale = 1.0/getDeltaT();
   for (int i = 0; i < numberOfAtoms; ++i) {
       if (masses[i] != 0.0)
           for (int j = 0; j < 3; ++j) {
               velocities[i][j] = velocityScale*(xPrime[i][j] - atomCoordinates[i][j]);
               atomCoordinates[i][j] = xPrime[i][j];
           }
   }
   ReferenceVirtualSites::computePositions(system, atomCoordinates);
   incrementTimeStep();
}







