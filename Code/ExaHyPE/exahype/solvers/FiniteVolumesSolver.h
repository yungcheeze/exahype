/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released under the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 **/

#ifndef _EXAHYPE_SOLVERS_FINITE_VOLUMES_SOLVER_H_
#define _EXAHYPE_SOLVERS_FINITE_VOLUMES_SOLVER_H_

#include "exahype/solvers/Solver.h"

#include "exahype/records/FiniteVolumesCellDescription.h"

namespace exahype {
namespace solvers {
class FiniteVolumesSolver;
}  // namespace solvers
}  // namespace exahype

class exahype::solvers::FiniteVolumesSolver : public exahype::solvers::Solver {
 private:
  /**
   * Total number of unknowns in a cell.
   */
  int _unknownsPerCell;

  /**
   * Total number of unknowns per cell face.
   */
  int _unknownsPerFace;

  /**
   * Total number of unknowns per cell boundary.
   */
  int _unknownsPerCellBoundary;

  /**
   * Minimum time stamps of all patches.
   */
  double _minTimeStamp;

  /**
   * Minimum time step size of all patches.
   */
  double _minTimeStepSize;

  /**
   * Next minimum step size of all patches.
   */
  double _nextMinTimeStepSize;

#ifdef Parallel
  /**
   * Data messages per neighbour communication.
   * This information is required by the sendEmpty...(...)
   * method.
   */
  static const int DataMessagesPerNeighbourCommunication;
  /**
   * Data messages per fork/join communication.
   * This information is required by the sendEmpty...(...)
   * method.
   */
  static const int DataMessagesPerForkOrJoinCommunication;
  /**
   * Data messages per master worker communication.
   * This information is required by the sendEmpty...(...)
   * method.
   */
  static const int DataMessagesPerMasterWorkerCommunication;
#endif

 public:
  typedef exahype::DataHeap DataHeap;

  /**
   * Rank-local heap that stores FiniteVolumesCellDescription instances.
   *
   * \note This heap might be shared by multiple FiniteVolumesSolver instances
   * that differ in their solver number and other attributes.
   * @see solvers::Solver::RegisteredSolvers.
   */
  typedef peano::heap::PlainHeap<exahype::records::FiniteVolumesCellDescription>
      Heap;

  FiniteVolumesSolver(const std::string& identifier, int numberOfVariables,
                      int numberOfParameters, int nodesPerCoordinateAxis,
                      double maximumMeshSize,
                      exahype::solvers::Solver::TimeStepping timeStepping,
                      std::unique_ptr<profilers::Profiler> profiler =
                          std::unique_ptr<profilers::Profiler>(
                              new profilers::simple::NoOpProfiler("")));

  virtual ~FiniteVolumesSolver() {}

  // Disallow copy and assignment
  FiniteVolumesSolver(const FiniteVolumesSolver& other) = delete;
  FiniteVolumesSolver& operator=(const FiniteVolumesSolver& other) = delete;

  /**
   * @param luh is a pointer to 3^d pointers to doubles
   */
  virtual double stableTimeStepSize(
      double* luh[THREE_POWER_D],
      const tarch::la::Vector<DIMENSIONS, double>& dx) = 0;

  /**
   *
   */
  virtual void solutionAdjustment(
      double* luh, const tarch::la::Vector<DIMENSIONS, double>& center,
      const tarch::la::Vector<DIMENSIONS, double>& dx, double t, double dt) = 0;

  virtual bool hasToAdjustSolution(
      const tarch::la::Vector<DIMENSIONS, double>& center,
      const tarch::la::Vector<DIMENSIONS, double>& dx, double t) = 0;

  virtual exahype::solvers::Solver::RefinementControl refinementCriterion(
      const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center,
      const tarch::la::Vector<DIMENSIONS, double>& dx, double t,
      const int level) = 0;

  /**
   * @param luh is a pointer to 3^d pointers to doubles
   * @param dt Time step size that is to be used.
   * @param maxAdmissibleDt Maximum time step size that would have been
   *        possible. If maxAdmissibleDt<dt, then we know that no time
   *        step has been done.
   */
  virtual void solutionUpdate(double* luh[THREE_POWER_D],
                              const tarch::la::Vector<DIMENSIONS, double>& dx,
                              const double dt, double& maxAdmissibleDt) = 0;

  virtual double getMinTimeStamp() const override;

  /**
   * This operation returns the number of unknowns per cell located in
   * the interior of a cell.
   */
  int getUnknownsPerCell() const;

  /**
   * Run over all solvers and identify the minimal time step size.
   */
  virtual double getMinTimeStepSize() const override;

  virtual void updateNextTimeStepSize(double value) override;

  virtual void initInitialTimeStamp(double value) override;

  void synchroniseTimeStepping(const int cellDescriptionsIndex,
                               const int element) override;

  virtual void startNewTimeStep() override;

  virtual double getNextMinTimeStepSize() const override;

  virtual int tryGetElement(const int cellDescriptionsIndex,
                            const int solverNumber) const;

#ifdef Parallel
  /**
   * Sends all the cell descriptions at address \p
   * cellDescriptionsIndex to the rank \p toRank.
   */
  static void sendCellDescriptions(
      const int toRank, const int cellDescriptionsIndex,
      const peano::heap::MessageType& messageType,
      const tarch::la::Vector<DIMENSIONS, double>& x, const int level);

  /**
   * Sends an empty message to the rank \p toRank.
   */
  static void sendEmptyCellDescriptions(
      const int toRank, const peano::heap::MessageType& messageType,
      const tarch::la::Vector<DIMENSIONS, double>& x, const int level);

  /**
   * Receives cell descriptions from rank \p fromRank
   * and resets the data heap indices to -1.
   *
   * If a received cell description has the same
   * solver number as a cell description in the
   * array at address \p cellDescriptionsIndex,
   * we merge the metadata (time stamps, time step size)
   * of both cell descriptions.
   *
   * If no cell description in the array at address
   * \p cellDescriptionsIndex can be found with the
   * same solver number than a received cell description,
   * we push the received cell description to
   * the back of the array at address \p cellDescriptions
   * Index.
   *
   * This operation is intended to be used in combination
   * with the solver method mergeWithWorkerOrMasterDataDueToForkOrJoin(...).
   * Here, we would merge first the cell descriptions sent by the master and
   * worker
   * and then merge the data that is sent out right after.
   */
  static void mergeCellDescriptionsWithRemoteData(
      const int fromRank, const int cellDescriptionsIndex,
      const peano::heap::MessageType& messageType,
      const tarch::la::Vector<DIMENSIONS, double>& x, const int level);

  /**
   * Drop cell descriptions received from \p fromRank.
   */
  static void dropCellDescriptions(
      const int fromRank, const peano::heap::MessageType& messageType,
      const tarch::la::Vector<DIMENSIONS, double>& x, const int level);

  /**
   * @deprecated
   */
  void sendToRank(int rank, int tag) override;

  /**
   * @deprecated
   */
  void receiveFromMasterRank(int rank, int tag) override;

  ///////////////////////////////////
  // NEIGHBOUR
  ///////////////////////////////////

  /**
   * Send solver data to neighbour rank.
   */
  void sendDataToNeighbour(const int toRank, const int cellDescriptionsIndex,
                           const int elementIndex,
                           const tarch::la::Vector<DIMENSIONS, int>& src,
                           const tarch::la::Vector<DIMENSIONS, int>& dest,
                           const tarch::la::Vector<DIMENSIONS, double>& x,
                           const int level) override {
    assertionMsg(false, "Please implement!");
  }

  /**
   * Send empty solver data to neighbour rank.
   */
  void sendEmptyDataToNeighbour(const int toRank,
                                const tarch::la::Vector<DIMENSIONS, int>& src,
                                const tarch::la::Vector<DIMENSIONS, int>& dest,
                                const tarch::la::Vector<DIMENSIONS, double>& x,
                                const int level) override;

  /**
   * Receive solver data from neighbour rank.
   *
   * \param[in] elementIndex Index of the ADERDGCellDescription
   *                         holding the data to send out in
   *                         the heap vector at \p cellDescriptionsIndex.
   */
  void mergeWithNeighbourData(const int fromRank,
                              const int cellDescriptionsIndex,
                              const int elementIndex,
                              const tarch::la::Vector<DIMENSIONS, int>& src,
                              const tarch::la::Vector<DIMENSIONS, int>& dest,
                              const tarch::la::Vector<DIMENSIONS, double>& x,
                              const int level) override {
    assertionMsg(false, "Please implement!");
  }

  /**
   * Drop solver data from neighbour rank.
   */
  void dropNeighbourData(const int fromRank,
                         const tarch::la::Vector<DIMENSIONS, int>& src,
                         const tarch::la::Vector<DIMENSIONS, int>& dest,
                         const tarch::la::Vector<DIMENSIONS, double>& x,
                         const int level) override;

  ///////////////////////////////////
  // FORK OR JOIN
  ///////////////////////////////////

  /**
   * Send solver data to master or worker rank. Read the data from
   * the cell description \p element in
   * the cell descriptions vector stored at \p
   * cellDescriptionsIndex.
   *
   * \param[in] element Index of the FiniteVolumesCellDescription
   *                    holding the data to send out in
   *                    the heap vector at \p cellDescriptionsIndex.
   */
  void sendDataToWorkerOrMasterDueToForkOrJoin(
      const int toRank, const int cellDescriptionsIndex, const int element,
      const tarch::la::Vector<DIMENSIONS, double>& x, const int level) override;

  /**
   * Send empty solver data to master or worker rank
   * due to fork or join.
   */
  void sendEmptyDataToWorkerOrMasterDueToForkOrJoin(
      const int toRank, const tarch::la::Vector<DIMENSIONS, double>& x,
      const int level) override;

  /**
   * Merge with solver data from master or worker rank
   * that was sent out due to a fork or join. Wrote the data to
   * the cell description \p element in
   * the cell descriptions vector stored at \p
   * cellDescriptionsIndex.
   */
  void mergeWithWorkerOrMasterDataDueToForkOrJoin(
      const int fromRank, const int cellDescriptionsIndex, const int element,
      const tarch::la::Vector<DIMENSIONS, double>& x, const int level) override;

  /**
   * Drop solver data from master or worker rank
   * that was sent out due to a fork or join.
   */
  void dropWorkerOrMasterDataDueToForkOrJoin(
      const int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
      const int level) override;

  ///////////////////////////////////
  // WORKER->MASTER
  ///////////////////////////////////

  /**
   * Send solver data to master rank. Read the data from
   * the cell description \p element in
   * the cell descriptions vector stored at \p
   * cellDescriptionsIndex.
   *
   * \param[in] element Index of the ADERDGCellDescription
   *                    holding the data to send out in
   *                    the heap vector at \p cellDescriptionsIndex.
   */
  void sendDataToMaster(const int masterRank, const int cellDescriptionsIndex,
                        const int element,
                        const tarch::la::Vector<DIMENSIONS, double>& x,
                        const int level) override;

  /**
   * Send empty solver data to master rank.
   */
  void sendEmptyDataToMaster(const int masterRank,
                             const tarch::la::Vector<DIMENSIONS, double>& x,
                             const int level) override;

  /**
   * Merge with solver data from worker rank.
   * Write the data to the cell description \p element in
   * the cell descriptions vector stored at \p
   * cellDescriptionsIndex.
   *
   * \param[in] element Index of the cell description
   *                    holding the data to send out in
   *                    the array with address \p cellDescriptionsIndex.
   */
  void mergeWithWorkerData(const int workerRank,
                           const int cellDescriptionsIndex, const int element,
                           const tarch::la::Vector<DIMENSIONS, double>& x,
                           const int level) override;

  /**
   * Drop solver data from worker rank.
   */
  void dropWorkerData(const int workerRank,
                      const tarch::la::Vector<DIMENSIONS, double>& x,
                      const int level) override;

  ///////////////////////////////////
  // MASTER->WORKER
  ///////////////////////////////////

  /**
   * Send solver data to worker rank. Read the data from
   * the cell description \p element in the cell descriptions
   * vector stored at \p cellDescriptionsIndex.
   *
   * \param[in] element Index of the ADERDGCellDescription
   *                    holding the data to send out in
   *                    the heap vector at \p cellDescriptionsIndex.
   */
  void sendDataToWorker(const int workerRank, const int cellDescriptionsIndex,
                        const int element,
                        const tarch::la::Vector<DIMENSIONS, double>& x,
                        const int level) override;

  /**
   * Send empty solver data to worker rank.
   */
  void sendEmptyDataToWorker(const int workerRank,
                             const tarch::la::Vector<DIMENSIONS, double>& x,
                             const int level) override;

  /**
   * Merge with solver data from master rank
   * that was sent out due to a fork or join. Write the data to
   * the cell description \p element in
   * the cell descriptions vector stored at \p
   * cellDescriptionsIndex.
   *
   * \param[in] element Index of the cell description
   *                    holding the data to send out in
   *                    the array with address \p cellDescriptionsIndex.
   */
  void mergeWithMasterData(const int masterRank,
                           const int cellDescriptionsIndex, const int element,
                           const tarch::la::Vector<DIMENSIONS, double>& x,
                           const int level) override;

  /**
   * Drop solver data from master rank.
   */
  void dropMasterData(const int masterRank,
                      const tarch::la::Vector<DIMENSIONS, double>& x,
                      const int level) override;
#endif

  /**
   * Returns a string representation of this solver.
   */
  virtual std::string toString() const;

  /**
   * Writes a string representation of this solver
   * to \p out.
   */
  virtual void toString(std::ostream& out) const;
};

#endif
