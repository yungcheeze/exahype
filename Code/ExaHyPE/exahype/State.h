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
 
#ifndef _EXAHYPE_STATE_H_
#define _EXAHYPE_STATE_H_

#include "exahype/records/State.h"
#include "peano/grid/State.h"

#include <vector>
#include <memory>

#include "peano/grid/Checkpoint.h"

namespace exahype {
class State;

/**
 * Forward declaration
 */
class Vertex;
/**
 * Forward declaration
 */
class Cell;

namespace repositories {
  /**
   * Forward declaration
   */
  class RepositoryArrayStack;
    class RepositorySTDStack;
  }
}



/**
 * Blueprint for solver state.
 *
 * This file has originally been created by the PDT and may be manually extended
 *to
 * the needs of your application. We do not recommend to remove anything!
 */
class exahype::State : public peano::grid::State<exahype::records::State> {
 private:
  typedef class peano::grid::State<exahype::records::State> Base;

  /**
   * Needed for checkpointing.
   */
  friend class exahype::repositories::RepositoryArrayStack;
  friend class exahype::repositories::RepositorySTDStack;

  void writeToCheckpoint(
      peano::grid::Checkpoint<Vertex, Cell>& checkpoint) const;
  void readFromCheckpoint(
      const peano::grid::Checkpoint<Vertex, Cell>& checkpoint);


  #ifdef Parallel
  /**
   * We need/use this field in the parallel mode, but we use it on the global
   * master only in operation updateRegularInitialGridRefinementStrategy(). We
   * memorise the idle ranks per lookup. If it has changed, we assume that some
   * load balancing is going on.
   */
  int _idleRanksAtLastLookup;
  #endif
 public:
  /**
   * Default Constructor
   *
   * This constructor is required by the framework's data container. Do not
   * remove it.
   */
  State();

  /**
   * Constructor
   *
   * This constructor is required by the framework's data container. Do not
   * remove it. It is kind of a copy constructor that converts an object which
   * comprises solely persistent attributes into a full attribute. This very
   * functionality is implemented within the super type, i.e. this constructor
   * has to invoke the correponsing super type's constructor and not the super
   * type standard constructor.
   */
  State(const Base::PersistentState& argument);

  /**
   * Merge this state with another state
   *
   * @todo Clarify which stuff has to be merged
   */
  void merge(const State& anotherState);
  ///@}

  /**
   * Becomes nop in the serial case.
   *
   * In the parallel case, we have three different grid refinement strategies
   * implemented in mappings::Refinement.
   *
   * - Default: Just refine the grid in vertexLastTime(). The touch last
   *   convention (we always refine in touchVertexLastTime; see Peano
   *   cookbook) ensure that the grid is built up iteration by iteration.
   * - Veto: We do not refine though the refinement criterion would ask us
   *   to do so. We set the veto four out of five traversals and thus delay
   *   the setup and allow the load balancing to keep pace.
   * - Aggressive: Once we see on the global master that no rank is idle
   *   anymore, we switch into aggressive and make all ranks build up the
   *   whole grid in one sweep.
   */
  void updateRegularInitialGridRefinementStrategy();

  /**
   * In the serial version of the code, this predicate always holds. In the
   * parallel case, it holds if and only if all ranks are already busy. As the
   * routine only may be used by the setup of the regular initial grid, it
   * thus is reasonable to invoke enforceRefine in the parallel case if the
   * result it true.
   *
   * If this operation returns refineInitialGridInCreationalEvents(), also
   * refineInitialGridInTouchVertexLastTime() should hold in the parallel
   * mode. Without MPI, the two always are different.
   *
   * Please consult the Peano cookbook (Sect. 6.3.2) for details/rationale.
   */
  bool refineInitialGridInCreationalEvents() const;

  /**
   *
   * Please consult the Peano cookbook (Sect. 6.3.2) for details/rationale.
   *
   * Means that the computational regular initial grid is to be refined in
   * touchVertexLastTime(), but it basically also means that you may refined
   * the grid though perhaps not in the creational routines.
   */
  bool refineInitialGridInTouchVertexLastTime() const;

  /**
   * TODO(Dominic): Add docu.
   */
  records::State::MergeMode getMergeMode() const {
    return _stateData.getMergeMode();
  }

  /**
   * TODO(Dominic): Add docu.
   */
  records::State::SendMode getSendMode() const {
    return _stateData.getSendMode();
  }

  /**
   * Merging and Sending contexts.
   * See both mappings for more details.
   */
  void switchToInitialConditionAndTimeStepSizeComputationContext() {
    switchToSolutionUpdateAndTimeStepSizeComputationContext();
  }

  void switchToPredictionAndTimeStepSizeComputationContext() {
    _stateData.setMergeMode(records::State::BroadcastAndMergeTimeStepData);
    _stateData.setSendMode (records::State::ReduceAndMergeTimeStepDataAndSendFaceData);
  }

  void switchToADERDGTimeStepContext() {
    _stateData.setFuseADERDGPhases(true);
    _stateData.setMergeMode(records::State::BroadcastAndMergeTimeStepDataAndMergeFaceData);
    _stateData.setSendMode (records::State::ReduceAndMergeTimeStepDataAndSendFaceData);
  }

  void switchToPredictionRerunContext() {
    switchToPredictionContext();
  }

  void switchToNeighbourDataMergingContext() {
    _stateData.setMergeMode(records::State::MergeFaceData);
    _stateData.setSendMode (records::State::SendNothing);
  }

  void switchToSolutionUpdateAndTimeStepSizeComputationContext() {
    _stateData.setFuseADERDGPhases(false);
    _stateData.setMergeMode(records::State::MergeNothing);
    _stateData.setSendMode (records::State::ReduceAndMergeTimeStepData);
  }

  void switchToPredictionContext() {
    _stateData.setFuseADERDGPhases(false);
    _stateData.setMergeMode(records::State::BroadcastAndMergeTimeStepData);
    _stateData.setSendMode (records::State::SendFaceData);
  }

  void setStabilityConditionOfOneSolverWasViolated(bool state) {
    _stateData.setStabilityConditionOfOneSolverWasViolated(state);
  }

  bool stabilityConditionOfOneSolverWasViolated() {
    return _stateData.getStabilityConditionOfOneSolverWasViolated();
  }

  void setFuseADERDGPhases(bool state) {
    _stateData.setFuseADERDGPhases(state);
  }

  bool fuseADERDGPhases() {
    return _stateData.getFuseADERDGPhases();
  }

  void setTimeStepSizeWeightForPredictionRerun(double value) {
    _stateData.setTimeStepSizeWeightForPredictionRerun(value);
  }

  double getTimeStepSizeWeightForPredictionRerun() {
    return _stateData.getTimeStepSizeWeightForPredictionRerun();
  }
};

#endif
