/*
 * ReferenceIndirectReconstructionDynamics.h
 *
 *  Created on: 2 Nov 2020
 *      Author: hannesvdc
 */

#ifndef PLATFORMS_REFERENCE_INCLUDE_REFERENCEINDIRECTRECONSTRUCTIONDYNAMICS_H_
#define PLATFORMS_REFERENCE_INCLUDE_REFERENCEINDIRECTRECONSTRUCTIONDYNAMICS_H_


#include "ReferenceDynamics.h"
#include "openmm/ReactionCoordinate.h"

namespace OpenMM {

class ReferenceIndirectReconstructionDynamics : public ReferenceDynamics {

   private:

      std::vector<OpenMM::Vec3> xPrime, macroVariable;
      double lambda;
      ReactionCoordinate* reactionCoordinate;

   public:

      /**---------------------------------------------------------------------------------------

         Constructor

         @param numberOfAtoms  number of atoms
         @param deltaT         delta t for dynamics
         @param temperature    temperature
         @param lambda         strength of the biasing potential

         --------------------------------------------------------------------------------------- */

      ReferenceIndirectReconstructionDynamics(int numberOfAtoms, double deltaT, double temperature, ReactionCoordinate* rc,  double lambda);

      /**---------------------------------------------------------------------------------------

         Destructor

         --------------------------------------------------------------------------------------- */

       ~ReferenceIndirectReconstructionDynamics();

      /**---------------------------------------------------------------------------------------

         Update

         @param system              the System to be integrated
         @param atomCoordinates     atom coordinates
         @param velocities          velocities
         @param forces              forces
         @param masses              atom masses
         @param tolerance           the constraint tolerance

         --------------------------------------------------------------------------------------- */

      void update(const OpenMM::System& system, std::vector<OpenMM::Vec3>& atomCoordinates,
                  std::vector<OpenMM::Vec3>& velocities, std::vector<OpenMM::Vec3>& forces, std::vector<double>& masses, double tolerance);

      void setMacroscopicVariable(std::vector<OpenMM::Vec3> z) {
    	  	  macroVariable = z;
      }

      std::vector<OpenMM::Vec3> getMacroscopicVariable() const {
    	  	  return macroVariable;
      }

      double getLambda() const {
    	      return lambda;
      }
};

} // namespace OpenMM


#endif /* PLATFORMS_REFERENCE_INCLUDE_REFERENCEINDIRECTRECONSTRUCTIONDYNAMICS_H_ */
