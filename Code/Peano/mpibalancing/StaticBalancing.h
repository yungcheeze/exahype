// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www.peano-framework.org
#ifndef _MPI_BALANCING_STATIC_BALANCING_H_
#define _MPI_BALANCING_STATIC_BALANCING_H_

#include "peano/parallel/loadbalancing/OracleForOnePhase.h"
#include "tarch/logging/Log.h"


#include <map>
#include <set>


namespace mpibalancing {
  class StaticBalancing;
}


/**
 * Static Balancing
 *
 * This balancing identifies the critical path locally, i.e. if the oracle is
 * told by its master that it shall fork at least once, it identifies its worker
 * with the biggest load and makes this one fork as well. Please note that the
 * critical path can be multiple ranks that are forked at the same time.
 *
 * !!! getCommandForWorker()
 *
 * getCommandForWorker() here is pretty simple and evaluates basically
 * _forksToConduct and _criticalWorker identifying whether a worker shall
 * fork (only if it is member of the set _criticalWorker) and how often it
 * shall fork (_forksToConduct). There are few special cases besides this:
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
 * This operation triggers the critical path analysis which consists of two
 * steps:
 *
 * - Determine the critial path
 * - Compute the number of forks along the critical paths
 *
 * If its command from its master is unequal to fork, there's no analysis to
 * be done, i.e. we evaluate the critical path recursively which works as we
 * are working in a tree environment. Otherwise, the two operations map onto
 * the routines identifyCriticalPathes() and
 * computeMaxForksOnCriticalWorker().
 *
 * If we shall fork, we first identify the biggest and the smallest local
 * workload. The local workload is always greater or equal to 1.
 * If the worker (the one directly below the global
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
 *
 * @image html StaticBalancing.png
 * @author Tobias Weinzierl
 */
class mpibalancing::StaticBalancing: public peano::parallel::loadbalancing::OracleForOnePhase {
  public:
    /**
     * Constructor
     *
     * @param joinsAllowed The static balancing does, as the name suggests, not
     *                     join any partitions throughout the computation. Only
     *                     if grid erases do not pass through due to the domain
     *                     decomposition, it supports joins if this flag is set
     *                     true.
     *
     * @param coarsestRegularInnerAndOuterGridLevel If you use this balancing,
     *                     the grid is refined regularly up to level
     *                     coarsestRegularInnerAndOuterGridLevel - indepenent
     *                     of whether grid elements are inside or outside of
     *                     the domain. Too regular grids facilitate a
     *                     proper load balancing in several cases.
     *
     * @param finestAdministrativeLevel Up to this level, all reanks try to yield
     *                     all workers regardless of the grammar. Use 0 to switch
     *                     off the total worker deployment everywhere besides for
     *                     the root node and its first worker (here, the total
     *                     deployment is required to make the overall balancing
     *                     work).
     */
    StaticBalancing( bool joinsAllowed, int coarsestRegularInnerAndOuterGridLevel = 3 );

    virtual ~StaticBalancing();

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
     * Runs an analysis on the _weightMap. This operation returns consistent
     * data if it is used within receivedStartCommand() and if all children
     * do restrict. It then holds all the weights from the previous iteration.
     * Otherwise, some results might be older (two or three iterations ago,
     * e.g.) or some values might also have been updated in the current
     * traversal.
     */
    double getMaximumWeightOfWorkers() const;

    /**
     * @see getMaximumWeightOfWorkers()
     * @see Class documentation clarifying why we have to exclude weights
     *      smaller than one explicitly.
     */
    double getMinimumWeightOfWorkers() const;

    /**
     * Run through all the workers and insert those that have a critical weight
     * into the _criticalWorker set. The latter is cleared upon entering the
     * routine.
     *
     * If the command received from the master does not allow us to fork
     * further, the set of critical ranks remains empty upon termination.
     *
     * If the command from the master is to become an administrative rank, we
     * insert all workers as critical workers.
     *
     * If there is no real minimum weight, we assume that we have run into a
     * local minimum and thus insert all workers into the fork list. See remark
     * in the class documentation.
     *
     * @pre _forksToConduct has to be set to command received from master
     */
    void identifyCriticalPathes( int commandFromMaster );

    void computeMaxForksOnCriticalWorker( const int commandFromMaster );

    /**
     * Logging device
     */
    static tarch::logging::Log  _log;

    /**
     * If a fork failed, all the oracles should stop to ask for further forks.
     * Wouldn't make sense and just slow down the application.
     */
    static bool                 _forkHasFailed;

    const int                   _coarsestRegularInnerAndOuterGridLevel;

    /**
     * Global flag set at construction time.
     */
    const bool                  _joinsAllowed;
    
    /**
     * Holds for each worker its local weight.
     */
    std::map<int,double>        _weightMap;

    /**
     * Bookkeeps for each worker whether this one couldn't erase due to the
     * domain decomposition.
     */
    std::map<int,bool>          _workerCouldNotEraseDueToDecomposition;

    /**
     * Set of critical workers
     *
     * This set might be empty, but it can also contain multiple critical
     * workers. All critical workers will receive fork commands if Peano
     * asks for a load balancing command next. What type of fork command
     * they receive is held in _forksToConduct.
     */
    std::set<int>               _criticalWorker;

    /**
     * Determines the number of forks for a worker along the critical path.
     */
    int                         _maxForksOnCriticalWorker;
};



#endif
