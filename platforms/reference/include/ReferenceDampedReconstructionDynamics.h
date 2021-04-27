/*
 * ReferenceIndirectReconstructionDynamics.h
 *
 *  Created on: 2 Nov 2020
 *      Author: hannesvdc
 */

#ifndef PLATFORMS_REFERENCE_INCLUDE_REFERENCEDAMPEDRECONSTRUCTIONDYNAMICS_H_
#define PLATFORMS_REFERENCE_INCLUDE_REFERENCEDAMPEDRECONSTRUCTIONDYNAMICS_H_


#include "ReferenceDynamics.h"
#include "openmm/ReactionCoordinate.h"

namespace OpenMM {

class ReferenceDampedReconstructionDynamics : public ReferenceDynamics {

   private:

      std::vector<OpenMM::Vec3> xPrime, macroVariable;
      double lambda;
      double gamma;
      ReactionCoordinate* reactionCoordinate;

   public:

      /**---------------------------------------------------------------------------------------

         Constructor

         @param numberOfAtoms  number of atoms
         @param deltaT         delta t for dynamics
         @param temperature    temperature
         @param lambda         strength of the biasing potential
         @param gamma           the damping coefficient

         --------------------------------------------------------------------------------------- */

    ReferenceDampedReconstructionDynamics(int numberOfAtoms, double deltaT, double temperature, ReactionCoordinate* rc,  double lambda, double gamma);

      /**---------------------------------------------------------------------------------------

         Destructor

         --------------------------------------------------------------------------------------- */

       ~ReferenceDampedReconstructionDynamics();

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
    
    double getGamma() const {
        return gamma;
    }
};

} // namespace OpenMM


#endif /* PLATFORMS_REFERENCE_INCLUDE_REFERENCEINDIRECTRECONSTRUCTIONDYNAMICS_H_ */
