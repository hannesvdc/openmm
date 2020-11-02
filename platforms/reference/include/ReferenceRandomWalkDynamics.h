/*
 * ReferenceRandomWalkDynamics.h
 *
 *  Created on: 28 Oct 2020
 *      Author: hannesvdc
 */

#ifndef PLATFORMS_REFERENCE_INCLUDE_REFERENCERANDOMWALKDYNAMICS_H_
#define PLATFORMS_REFERENCE_INCLUDE_REFERENCERANDOMWALKDYNAMICS_H_


#include "ReferenceDynamics.h"

namespace OpenMM {

class ReferenceRandomWalkDynamics : public ReferenceDynamics {

   private:

      std::vector<OpenMM::Vec3> xPrime;

   public:

      /**---------------------------------------------------------------------------------------

         Constructor

         @param numberOfAtoms  number of atoms
         @param deltaT         delta t for dynamics
         @param friction       friction coefficient
         @param temperature    temperature

         --------------------------------------------------------------------------------------- */

       ReferenceRandomWalkDynamics(int numberOfAtoms, double deltaT, double temperature);

      /**---------------------------------------------------------------------------------------

         Destructor

         --------------------------------------------------------------------------------------- */

       ~ReferenceRandomWalkDynamics();

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

};

} // namespace OpenMM


#endif /* PLATFORMS_REFERENCE_INCLUDE_REFERENCERANDOMWALKDYNAMICS_H_ */
