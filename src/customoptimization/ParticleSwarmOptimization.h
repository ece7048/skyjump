#ifndef PARTICLE_SWARM_OPTIMIZER_H
#define PARTICLE_SWARM_OPTIMIZER_H

#include <SimTKcommon.h>

#include <simmath/internal/common.h>
#include <simmath/Optimizer.h>
#include <simmath/internal/OptimizerRep.h>

namespace SimTK
{
    class ParticleSwarmOptimizer : public Optimizer::OptimizerRep
    {
    public:
        ~ParticleSwarmOptimizer()
        {
            delete [] g_U;
            delete [] g_L;
            delete [] mult_x_L;
            delete [] mult_x_U;
            delete [] mult_g;
        }

        ParticleSwarmOptimizer(const OptimizerSystem& sys);

        Real optimize(Vector &results) override;
        OptimizerRep* clone() const override;

        OptimizerAlgorithm getAlgorithm() const override
        {
            return InteriorPoint;
        }

    private:
        Real         *mult_x_L;
        Real         *mult_x_U;
        Real         *mult_g;
        Real         *g_L;
        Real         *g_U;
        bool          firstOptimization;
    };
} // namespace SimTK

#endif 