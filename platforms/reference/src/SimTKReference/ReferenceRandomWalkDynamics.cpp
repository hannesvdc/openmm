#include <cstring>
#include <sstream>

#include "SimTKOpenMMUtilities.h"
#include "ReferenceRandomWalkDynamics.h"
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

ReferenceRandomWalkDynamics::ReferenceRandomWalkDynamics(int numberOfAtoms,
                                                          double deltaT,
                                                          double temperature) :
           ReferenceDynamics(numberOfAtoms, deltaT, temperature) {

   xPrime.resize(numberOfAtoms);
}

/**---------------------------------------------------------------------------------------

   ReferenceBrownianDynamics destructor

   --------------------------------------------------------------------------------------- */

ReferenceRandomWalkDynamics::~ReferenceRandomWalkDynamics() {
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

void ReferenceRandomWalkDynamics::update(const OpenMM::System& system, vector<Vec3>& atomCoordinates,
                                          vector<Vec3>& velocities,
                                          vector<Vec3>& forces, vector<double>& masses, double tolerance) {

   // Perform the integration.
   int numberOfAtoms = system.getNumParticles();
   const double noiseAmplitude = sqrt(2.0*BOLTZ*getTemperature()*getDeltaT());
   for (int i = 0; i < numberOfAtoms; ++i) {
       if (masses[i] != 0.0)
           for (int j = 0; j < 3; ++j) {
               xPrime[i][j] = atomCoordinates[i][j] + noiseAmplitude*SimTKOpenMMUtilities::getNormallyDistributedRandomNumber();
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




