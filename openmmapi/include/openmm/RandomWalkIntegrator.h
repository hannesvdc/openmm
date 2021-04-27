/*
 * RandomWalkIntegrator.h
 *
 *  Created on: 28 Oct 2020
 *      Author: hannesvdc
 */

#ifndef OPENMM_RANDOMWALKINTEGRATOR_H_
#define OPENMM_RANDOMWALKINTEGRATOR_H_

#include "Integrator.h"
#include "openmm/Kernel.h"
#include "internal/windowsExport.h"

namespace OpenMM {

/**
 * This is an Integrator which simulates a System using Brownian dynamics.
 */
class OPENMM_EXPORT RandomWalkIntegrator : public Integrator {
public:
	/**
	 * Create a RandomWalkIntegrator.
	 *
	 * @param temperature    the temperature of the heat bath (in Kelvin)
	 * @param stepSize       the step size with which to integrate the system (in picoseconds)
	 */
	RandomWalkIntegrator(double temperature, double stepSize, double period=-1.);
	/**
	 * Get the temperature of the heat bath (in Kelvin).
	 *
	 * @return the temperature of the heat bath (in Kelvin).
	 */
	double getTemperature() const {
		return temperature;
	}
	/**
	 * Set the temperature of the heat bath (in Kelvin).
	 *
	 * @param temp    the temperature of the heat bath, measured in Kelvin.
	 */
	void setTemperature(double temp) {
		temperature = temp;
	}
	/**
	 * Get the random number seed.  See setRandomNumberSeed() for details.
	 */
	int getRandomNumberSeed() const {
		return randomNumberSeed;
	}
	/**
	 * Set the random number seed.  The precise meaning of this parameter is undefined, and is left up
	 * to each Platform to interpret in an appropriate way.  It is guaranteed that if two simulations
	 * are run with different random number seeds, the sequence of random forces will be different.  On
	 * the other hand, no guarantees are made about the behavior of simulations that use the same seed.
	 * In particular, Platforms are permitted to use non-deterministic algorithms which produce different
	 * results on successive runs, even if those runs were initialized identically.
	 *
	 * If seed is set to 0 (which is the default value assigned), a unique seed is chosen when a Context
	 * is created from this Force. This is done to ensure that each Context receives unique random seeds
	 * without you needing to set them explicitly.
	 */
	void setRandomNumberSeed(int seed) {
		randomNumberSeed = seed;
	}

	/**
	 * Setup the MCMC sampler by computing the current state.
	 */
	void setupSampler();

	/**
	 * Pass whether the proposal move as accepted or not.
	 *
	 * @param acc     true if the proposal move was accepted, false if not
	 */
	void accepted(bool acc);

	/**
	 * Get the previous state.
	 */
//	State getPreviousState() const {
//		return prev_state;
//	}

	/**
	 * Advance a simulation through time by taking a series of time steps.
	 *
	 * @param steps   the number of time steps to take
	 */
	void step(int steps);
protected:
	/**
	 * This will be called by the Context when it is created.  It informs the Integrator
	 * of what context it will be integrating, and gives it a chance to do any necessary initialization.
	 * It will also get called again if the application calls reinitialize() on the Context.
	 */
	void initialize(ContextImpl& context);
	/**
	 * This will be called by the Context when it is destroyed to let the Integrator do any necessary
	 * cleanup.  It will also get called again if the application calls reinitialize() on the Context.
	 */
	void cleanup();
	/**
	 * Get the names of all Kernels used by this Integrator.
	 */
	std::vector<std::string> getKernelNames();
	/**
	 * Compute the kinetic energy of the system at the current time.
	 */
	double computeKineticEnergy();
	/**
	 * Computing kinetic energy for this integrator does not require forces.
	 */
	bool kineticEnergyRequiresForce() const;
private:
    double _period;
    double temperature;
    int randomNumberSeed;
    State prev_state;
    Kernel kernel;
};

} // namespace OpenMM



#endif /* OPENMM_RANDOMWALKINTEGRATOR_H_ */
