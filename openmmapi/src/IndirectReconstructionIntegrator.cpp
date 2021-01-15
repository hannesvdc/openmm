/*
 * IndirectReconstructionIntegrator.cpp
 *
 *  Created on: 2 Nov 2020
 *      Author: hannesvdc
 */


#include "openmm/IndirectReconstructionIntegrator.h"
#include "openmm/Context.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/kernels.h"
#include <string>
#include <iostream>

using namespace OpenMM;
using std::string;
using std::vector;

IndirectReconstructionIntegrator::IndirectReconstructionIntegrator(double temperature, double lambda, double stepSize, ReactionCoordinate* rc) {
    setTemperature(temperature);
    setStepSize(stepSize);
    setLambda(lambda);
    setConstraintTolerance(1e-5);
    setRandomNumberSeed(0);
    setReactionCoordinate(rc);
}

void IndirectReconstructionIntegrator::initialize(ContextImpl& contextRef) {
	if (owner != NULL && &contextRef.getOwner() != owner)
		throw OpenMMException("This Integrator is already bound to a context");
	context = &contextRef;
	owner = &contextRef.getOwner();
	kernel = context->getPlatform().createKernel(IntegrateIndirectReconstructionStepKernel::Name(), contextRef);
	kernel.getAs<IntegrateIndirectReconstructionStepKernel>().initialize(contextRef.getSystem(), *this);
}

void IndirectReconstructionIntegrator::cleanup() {
	kernel = Kernel();
}

vector<string> IndirectReconstructionIntegrator::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(IntegrateIndirectReconstructionStepKernel::Name());
    return names;
}

double IndirectReconstructionIntegrator::computeKineticEnergy() {
    return kernel.getAs<IntegrateIndirectReconstructionStepKernel>().computeKineticEnergy(*context, *this);
}

bool IndirectReconstructionIntegrator::kineticEnergyRequiresForce() const {
    return false;
}

void IndirectReconstructionIntegrator::setupSampler() {
	context->updateContextState();
	context->calcForcesAndEnergy(true, true, getIntegrationForceGroups());
	prev_state = owner->getState(State::Positions | State::Forces | State::Energy);
}

void IndirectReconstructionIntegrator::accepted(bool acc) {
	if (acc) {
		prev_state = owner->getState(State::Positions | State::Forces | State::Energy);
	} else {
		owner->setState(prev_state);
		context->updateContextState();
	}
}

void IndirectReconstructionIntegrator::step(int steps) {
	if (context == NULL)
	    throw OpenMMException("This Integrator is not bound to a context!");

	context->updateContextState();
	for (int i = 0; i < steps; ++i) {
	    kernel.getAs<IntegrateIndirectReconstructionStepKernel>().execute(*context, *this);
	    context->updateContextState();
	    context->calcForcesAndEnergy(true, true, getIntegrationForceGroups());
	}
}
