/*
 * IndirectReconstructionIntegrator.cpp
 *
 *  Created on: 2 Nov 2020
 *      Author: hannesvdc
 */


#include "openmm/DampedReconstructionIntegrator.h"
#include "openmm/Context.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/kernels.h"
#include <string>
#include <iostream>
#include <time.h>

using namespace OpenMM;
using std::string;
using std::vector;


DampedReconstructionIntegrator::DampedReconstructionIntegrator(double temperature, double lambda, double stepSize, double gamma, ReactionCoordinate* rc) {
    setTemperature(temperature);
    setStepSize(stepSize);
    setLambda(lambda);
    setGamma(gamma);
    setConstraintTolerance(1e-5);
    setRandomNumberSeed(0);
    setReactionCoordinate(rc);
    beta =  1./(8.3145*temperature/1000.);
}

void DampedReconstructionIntegrator::initialize(ContextImpl& contextRef) {
	if (owner != NULL && &contextRef.getOwner() != owner)
		throw OpenMMException("This Integrator is already bound to a context");
	context = &contextRef;
	owner = &contextRef.getOwner();
	kernel = context->getPlatform().createKernel(IntegrateDampedReconstructionStepKernel::Name(), contextRef);
	kernel.getAs<IntegrateDampedReconstructionStepKernel>().initialize(contextRef.getSystem(), *this);
    kernel.getAs<IntegrateDampedReconstructionStepKernel>().setLambda(lambda);
    kernel.getAs<IntegrateDampedReconstructionStepKernel>().setGamma(gamma);
    kernel.getAs<IntegrateDampedReconstructionStepKernel>().setReactionCoordinate(reactionCoordinate);
}

void DampedReconstructionIntegrator::cleanup() {
	kernel = Kernel();
}

vector<string> DampedReconstructionIntegrator::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(IntegrateDampedReconstructionStepKernel::Name());
    return names;
}

double DampedReconstructionIntegrator::computeKineticEnergy() {
    return kernel.getAs<IntegrateDampedReconstructionStepKernel>().computeKineticEnergy(*context, *this);
}

bool DampedReconstructionIntegrator::kineticEnergyRequiresForce() const {
    return false;
}

void DampedReconstructionIntegrator::setMacroscopicVariable(std::vector<Vec3> z) {
    macroVariable = z;
    kernel.getAs<IntegrateDampedReconstructionStepKernel>().setMacroscopicVariable(z);
}

void DampedReconstructionIntegrator::setupSampler() {
	context->updateContextState();
	context->calcForcesAndEnergy(true, true, getIntegrationForceGroups());
	prev_state = owner->getState(State::Positions | State::Forces | State::Energy);
    std::srand( (unsigned)time( NULL ) );
}

void DampedReconstructionIntegrator::accepted(bool acc) {
	if (acc) {
		prev_state = owner->getState(State::Positions | State::Forces | State::Energy);
	} else {
		owner->setState(prev_state);
		context->updateContextState();
	}
}

void DampedReconstructionIntegrator::step(int steps) {
	if (context == NULL)
	    throw OpenMMException("This Integrator is not bound to a context!");

	context->updateContextState();
    context->calcForcesAndEnergy(true, true, getIntegrationForceGroups());
	for (int i = 0; i < steps; ++i) {
	    kernel.getAs<IntegrateDampedReconstructionStepKernel>().execute(*context, *this);
        context->updateContextState();
	    context->calcForcesAndEnergy(true, true, getIntegrationForceGroups());
    }
}
