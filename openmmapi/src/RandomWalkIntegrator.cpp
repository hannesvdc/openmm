/*
 * RandomWalkIntegrator.cpp
 *
 *  Created on: 28 Oct 2020
 *      Author: hannesvdc
 */

#include "openmm/RandomWalkIntegrator.h"
#include "openmm/Context.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/kernels.h"
#include <string>
#include <iostream>

using namespace OpenMM;
using std::string;
using std::vector;

RandomWalkIntegrator::RandomWalkIntegrator(double temperature, double stepSize, double period) {
    setTemperature(temperature);
    setStepSize(stepSize);
    setConstraintTolerance(1e-5);
    setRandomNumberSeed(0);
    _period = period;
}

void RandomWalkIntegrator::initialize(ContextImpl& contextRef) {
	if (owner != NULL && &contextRef.getOwner() != owner)
		throw OpenMMException("This Integrator is already bound to a context");
	context = &contextRef;
	owner = &contextRef.getOwner();
	kernel = context->getPlatform().createKernel(IntegrateRandomWalkStepKernel::Name(), contextRef);
	kernel.getAs<IntegrateRandomWalkStepKernel>().initialize(contextRef.getSystem(), *this);
    kernel.getAs<IntegrateRandomWalkStepKernel>().setPeriod(_period);
}

void RandomWalkIntegrator::cleanup() {
	kernel = Kernel();
}

vector<string> RandomWalkIntegrator::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(IntegrateRandomWalkStepKernel::Name());
    return names;
}

double RandomWalkIntegrator::computeKineticEnergy() {
    return kernel.getAs<IntegrateRandomWalkStepKernel>().computeKineticEnergy(*context, *this);
}

bool RandomWalkIntegrator::kineticEnergyRequiresForce() const {
    return false;
}

void RandomWalkIntegrator::setupSampler() {
	context->updateContextState();
	double energy = context->calcForcesAndEnergy(true, true, getIntegrationForceGroups());
	prev_state = owner->getState(State::Positions | State::Energy);
}

void RandomWalkIntegrator::accepted(bool acc) {
	if (acc) {
		prev_state = owner->getState(State::Positions | State::Energy);
	} else {
		owner->setState(prev_state);
		context->updateContextState();
	}
}

void RandomWalkIntegrator::step(int steps) {
	if (context == NULL)
	    throw OpenMMException("This Integrator is not bound to a context!");

	for (int i = 0; i < steps; ++i) {
	    context->updateContextState();
	    kernel.getAs<IntegrateRandomWalkStepKernel>().execute(*context, *this);
	}

	context->updateContextState();
    context->calcForcesAndEnergy(false, true, getIntegrationForceGroups());
}



