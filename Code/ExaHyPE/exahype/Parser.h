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
 
#ifndef EXAHYPE_PARSER
#define EXAHYPE_PARSER

namespace exahype {
class Parser;
}

#include <iostream>

#include <vector>
#include <map>

#include "peano/utils/Globals.h"
#include "tarch/la/Vector.h"
#include "tarch/logging/Log.h"

#include "exahype/solvers/Solver.h"

/**
 * ExaHyPE command line parser
 *
 * @author Tobias Weinzierl
 */
class exahype::Parser {
 private:
  static tarch::logging::Log _log;

  std::vector<std::string> _tokenStream;

  /*
   * Helper map for converting strings to types.
   */
  std::map<std::string,exahype::solvers::Solver::Type> _identifier2Type;

  /*
   * Helper map for converting strings to types.
   */
  std::map<std::string,exahype::solvers::Solver::TimeStepping> _identifier2TimeStepping;

  /**
   * \return "notoken" if not found.
   */
  std::string getTokenAfter(std::string token,
                            int additionalTokensToSkip = 0) const;
  std::string getTokenAfter(std::string token0, std::string token1,
                            int additionalTokensToSkip = 0) const;
  std::string getTokenAfter(std::string token0, int occurance0,
                            int additionalTokensToSkip) const;
  std::string getTokenAfter(std::string token0, int occurance0,
                            std::string token1, int occurance1,
                            int additionalTokensToSkip = 0) const;
 public:
  Parser();
  virtual ~Parser() {}

  // Disallow copy and assignment
  Parser(const Parser& other) = delete;
  Parser& operator=(const Parser& other) = delete;

  enum MulticoreOracleType {
    Dummy,
    Autotuning,
    GrainSizeSampling
    // evtl. spaeter mal InvadeSHM
  };

  void readFile(const std::string& filename);

  bool isValid() const;

  /**
   * \return How many threads is the code supposed to use?
   */
  int getNumberOfThreads() const;

  tarch::la::Vector<DIMENSIONS, double> getDomainSize() const;
  tarch::la::Vector<DIMENSIONS, double> getBoundingBoxSize() const;

  tarch::la::Vector<DIMENSIONS, double> getOffset() const;

  std::string getMulticorePropertiesFile() const;

  MulticoreOracleType getMulticoreOracleType() const;

  double getSimulationEndTime() const;

  /**
   * \return Indicates if the user has chosen the fused ADER-DG time stepping variant.
   */
  bool getFuseAlgorithmicSteps() const;

  /**
   * \return Time step size underestimation factor for the fused ADER-DG time stepping variant.
   */
  double getFuseAlgorithmicStepsFactor() const;

  /**
   * \return The type of a solver.
   */
  exahype::solvers::Solver::Type getType(int solverNumber) const;

  /**
   * \return The identifier of a solver.
   */
  std::string getIdentifier(int solverNumber) const;

  /**
   * \return The number of state variables of a solver.
   */
  int getVariables(int solverNumber) const;

  /**
   * \return The number of parameters of a solver, e.g. material values etc.
   */
  int getParameters(int solverNumber) const;

  /**
   * \return The order of the ansatz polynomials of a solver.
   */
  int getOrder(int solverNumber) const;

  /**
   * \return The maximum extent in each coordinate direction a cell is allowed to have.
   */
  tarch::la::Vector<DIMENSIONS, double> getMaximumMeshSize(int solverNumber) const;

  /**
   * Prints a summary of the parameters read in for a solver.
   */
  void logSolverDetails(int solverNumber) const;

  /**
   * \return The time stepping mode of a solver.
   */
  exahype::solvers::Solver::TimeStepping getTimeStepping(int solverNumber) const;

  double getFirstSnapshotTimeForPlotter(int solverNumber,
                                        int plotterNumber) const;
  double getRepeatTimeForPlotter(int solverNumber, int plotterNumber) const;
  std::string getSelectorForPlotter(int solverNumber, int plotterNumber) const;
  std::string getIdentifierForPlotter(int solverNumber,
                                      int plotterNumber) const;
  std::string getFilenameForPlotter(int solverNumber, int plotterNumber) const;

  std::string getProfilerIdentifier() const;
  std::string getMetricsIdentifierList() const;
};

#endif
