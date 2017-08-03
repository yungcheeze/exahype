// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www.peano-framework.org
#ifndef _MPIBALANCING_SFC_DIFFUSION_NODE_POOL_STRATEGY_TEST_H_
#define _MPIBALANCING_SFC_DIFFUSION_NODE_POOL_STRATEGY_TEST_H_


#include "tarch/tests/TestCase.h"
#include "peano/grid/tests/records/TestCell.h"
#include "peano/grid/Cell.h"


#include "mpibalancing/SFCDiffusionNodePoolStrategy.h"


namespace mpibalancing {
  namespace tests {
    class SFCDiffusionNodePoolStrategyTest;
  }
}


/**
 *
 * @author Tobias Weinzierl
 */
class mpibalancing::tests::SFCDiffusionNodePoolStrategyTest: public tarch::tests::TestCase {
  private:
    static tarch::logging::Log _log;

    SFCDiffusionNodePoolStrategy createSetupWith4Nodes() const;

    /**
     * @see testPrimaryNodeDeliveryWith4Nodes() for sequence of requests
     */
    SFCDiffusionNodePoolStrategy::RequestQueue createQueueTriggeredByWorkersOfFirstWorkerWith4Nodes() const;

    void testFirstWorkerDelivery();
    void testPrimaryNodeDeliveryWith4Nodes();
    void testSecondaryNodeDeliveryWith4Nodes();

  public:
    SFCDiffusionNodePoolStrategyTest();

    virtual ~SFCDiffusionNodePoolStrategyTest();

    virtual void run();
};

#endif
