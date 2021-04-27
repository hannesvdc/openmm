/*
 * ReactionCoordinate.h
 *
 *  Created on: 5 Nov 2020
 *      Author: hannesvdc
 */

#ifndef OPENMM_REACTIONCOORDINATE_H_
#define OPENMM_REACTIONCOORDINATE_H_

#include <vector>
#include "Vec3.h"

#include "internal/windowsExport.h"

namespace OpenMM {

class OPENMM_EXPORT ReactionCoordinate {
public:
	ReactionCoordinate();

	virtual ~ReactionCoordinate();

    virtual std::vector<Vec3> value(const std::vector<Vec3>& x) = 0;

	virtual std::vector<Vec3> gradMatMul(const std::vector<Vec3>& x, const std::vector<Vec3> &z) = 0;
    
    virtual double getBiasedEnergy(const std::vector<OpenMM::Vec3> &x,
                                   const std::vector<OpenMM::Vec3> &z) = 0;
};

} // namespace OpenMM

#endif /* OPENMM_REACTIONCOORDINATE_H_ */
