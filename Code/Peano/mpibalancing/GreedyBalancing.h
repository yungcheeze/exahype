// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www.peano-framework.org
#ifndef _MPI_BALANCING_GREEDY_BALANCING_H_
#define _MPI_BALANCING_GREEDY_BALANCING_H_

#include "peano/parallel/loadbalancing/OracleForOnePhase.h"
#include "tarch/logging/Log.h"


#include <map>
#include <set>


namespace mpibalancing {
  class GreedyBalancing;
}


/**
 * Greedy Balancing
 *
 * This balancing identifies the critical path locally, i.e. if the oracle is
 * told by its master that it shall fork at least once, it identifies its worker
 * with the biggest load and makes this one fork as well. Please note that the
 * critical path can be multiple ranks that are forked at the same time.
 *
 * This greedy strategy is slightly different to the greedy strategy provided
 * by Peano's kernel (the dummy/default implementation) where every rank tries
 * to fork as often as possible.
 *
 *
 * !!! getCommandForWorker()
 *
 * getCommandForWorker() here is pretty simple and evaluates basically
 * _forksToConduct and _criticalWorker identifying whether the worker shall
 * fork or not. There are few special cases besides this:
 *
 * - If the oracle is running on the global master and the weight table holds
 *   only one entry (i.e. the weight of the local rank), this is the very
 *   root of the logic mpi tree and calls its worker to fork all of its nodes,
 *   i.e. we make rank 0 and rank 1 fork all their children. This way, we can
 *   evaluate on rank 1 really which subtree is the biggest one. If we'd fork
 *   all subtrees besides one on rank 1, this would not be possible that easy,
 *   as it is complicated to identify the real local workload throughout a
 *   traversal.
 * - If the oracle is asked for the command for a worker who couldn't erase
 *   before due to the domain decomposition and if joins are globally allowed,
 *   we join this one to facilitate that erases pass through and the global
 *   grid becomes smaller.
 *
 *
 * !!! receivedStartCommand()
 *
 * This operation is the critical path analysis. If its command from above is
 * unequal to fork, there's no analysis to be done, i.e. we evaluate the
 * critical path recursively which works as we are working in a tree
 * environment.
 *
 * If we shall fork, we first identify the biggest and the smallest local
 * workload. For the latter, we have to be careful as the local workload is
 * by default set to 1. Now, if the worker (the one directly below the global
 * master, e.g.) forks all of its children, this local workload remains 1,
 * which is an unrealistic value. So we exclude that one explicitly from our
 * search.
 *
 * Once we have the minimum and maximum workload as well as the critical
 * worker, i.e. the one having the biggest workload (which might be the
 * local rank as well), we determine the number of forks to conduct. The
 * formula for this endeavour is as follows: We assume (at the beginning)
 * that k forks reduce the worker's workload by
 *
 * @f$ w \gets w \cdot \frac{3^d-k}{3^d} @f$
 *
 * This assumption is based upon a regular grid. If the critical weight shall
 * equal the smallest weight, we obtain
 *
 * @f$ w_{critical} \cdot \frac{3^d-k}{3^d} = w_{minimum} @f$
 * @f$ k = 3^d \left( 1-\frac{w_{minimum}}{w_{critical}} \right) @f$
 *
 * Finally, we ensure that the forks are bounded by @f$ 3^d-1 @f$ or reset
 * the critical worker if no forks are to be done.
 *
 * !!! Local Minima
 *
 * If we have p ranks, it can happen that the oracle achieves perfect balancing with
 * @f$ \hat p < p @f$ ranks already. If we have 50 ranks, 10 in 2d already are perfectly
 * balanced though then 40 ranks remain idle. In this case, the critical path analysis
 * runs into a local minimum. We identify such a minimum or we assume that we are in
 * such a minimum, if the critical weight differs from the minimum weight by less then
 * LocalMinimumThreshold percent.
 *
 * If we identify a local minimum, we set the numbers to fork manually to
 * @f$3^d-1 @f$ and label all workers as critical.
 * From now on, the formula from above that dermines fork numbers doesn't work anymore.
 * We cannot assume that any worker still has @f$ 3^d @f$ children to fork. As a
 * consequence, we scale the whole fork numbers down by a factor of 1.1. This factor is
 * empirically. Please note that the effective scaling is typically
 * @f$ 1.1 ^4 \approx 1.45 @f$ as the local minimum analysis will trigger several times
 * until the grid is balanced and can rebalance again. So whenever we encounter a local
 * minimum, we are more carefully with rebalancing.
 *
 * !!! Setup phase without grammar analysis
 *
 * In the setup phase, we do skip all the critical path analysis and just try to fork
 * all nodes and become a purely administrative rank. This is done for all nodes responsible
 * for cells above _finestAdministrativeLevel. The latter is set such that in a purely
 * regular case, it is basically the finest level where still all cells can be handled
 * by different ranks.
 *
 * @image html GreedyBalancing.png
 * @author Tobias Weinzierl
 */
class mpibalancing::GreedyBalancing: public peano::parallel::loadbalancing::OracleForOnePhase {
  public:
    /**
     * Constructor
     *
     * @param coarsestRegularInnerAndOuterGridLevel If you use this balancing,
     *                     the grid is refined regularly up to level
     *                     coarsestRegularInnerAndOuterGridLevel - indepenent
     *                     of whether grid elements are inside or outside of
     *                     the domain. Too regular grids facilitate a
     *                     proper load balancing in several cases.
     */
    GreedyBalancing( int coarsestLevelWithRealWork, int coarsestRegularInnerAndOuterGridLevel = 3 );

    virtual ~GreedyBalancing();

    /**
     * Analyse which workers are to be forked next.
     *
     * This operation is described in detail within the class documentation as
     * it realises the actual algorithm.
     */
    void receivedStartCommand(const int commandFromMaster ) override;

    int getCommandForWorker( int workerRank, bool forkIsAllowed, bool joinIsAllowed ) override;

    void receivedTerminateCommand(
      int     workerRank,
      double  workerNumberOfInnerVertices,
      double  workerNumberOfBoundaryVertices,
      double  workerNumberOfOuterVertices,
      double  workerNumberOfInnerCells,
      double  workerNumberOfOuterCells,
      int     workerMaxLevel,
      int     currentLevel,
      const tarch::la::Vector<DIMENSIONS,double>& boundingBoxOffset,
      const tarch::la::Vector<DIMENSIONS,double>& boundingBoxSize,
      bool    workerCouldNotEraseDueToDecomposition
    ) override;

    void plotStatistics() override;

    peano::parallel::loadbalancing::OracleForOnePhase* createNewOracle(int adapterNumber) const override;

    void forkFailed() override;
 
    int getCoarsestRegularInnerAndOuterGridLevel() const override;

  private:
    /**
     * Logging device
     */
    static tarch::logging::Log  _log;

    /**
     * If a fork failed, all the oracles should stop to ask for further forks.
     * Wouldn't make sense and just slow down the application.
     */
    static bool                 _forkHasFailed;

    const int                   _coarsestLevelWithRealWork;

    const int                   _coarsestRegularInnerAndOuterGridLevel;

    std::map<int,int>           _workersLevel;

    int                         _finestLevelToForkAggressively;
};



#endif
