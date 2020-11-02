#ifndef OPENMM_RANDOMWALK_INTEGRATOR_PROXY_H_
#define OPENMM_RANDOMWALK_INTEGRATOR_PROXY_H_

#include "openmm/serialization/XmlSerializer.h"

namespace OpenMM {

    class RandomWalkIntegratorProxy : public SerializationProxy {
    public:
     	RandomWalkIntegratorProxy();
        void serialize(const void* object, SerializationNode& node) const;
        void* deserialize(const SerializationNode& node) const;
    };

}

#endif /*OPENMM_RANDOMWALK_INTEGRATOR_PROXY_H_*/
