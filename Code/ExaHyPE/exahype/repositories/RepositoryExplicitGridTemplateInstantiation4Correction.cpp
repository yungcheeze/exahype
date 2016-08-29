#include "exahype/repositories/Repository.h"
#include "exahype/records/RepositoryState.h"

#include "exahype/State.h"
#include "exahype/Vertex.h"
#include "exahype/Cell.h"

#include "peano/grid/Grid.h"

#include "peano/stacks/CellArrayStack.h"
#include "peano/stacks/CellSTDStack.h"

#include "peano/stacks/VertexArrayStack.h"
#include "peano/stacks/VertexSTDStack.h"

 #include "exahype/adapters/AugmentedAMRGrid.h" 
 #include "exahype/adapters/PlotAugmentedAMRGrid.h" 
 #include "exahype/adapters/InitialConditionAndTimeStepSizeComputation.h" 
 #include "exahype/adapters/PredictionAndPlotAndTimeStepSizeComputation.h" 
 #include "exahype/adapters/PredictionAndTimeStepSizeComputation.h" 
 #include "exahype/adapters/GridErasing.h" 
 #include "exahype/adapters/ADERDGTimeStep.h" 
 #include "exahype/adapters/ADERDGTimeStepAndPlot.h" 
 #include "exahype/adapters/PredictionRerun.h" 
 #include "exahype/adapters/RiemannSolver.h" 
 #include "exahype/adapters/Prediction.h" 
 #include "exahype/adapters/Correction.h" 
 #include "exahype/adapters/CorrectionAndPlot.h" 
 #include "exahype/adapters/Plot.h" 


namespace peano {
  namespace grid {
    template class Grid<exahype::Vertex,exahype::Cell,exahype::State, peano::stacks::VertexArrayStack<exahype::Vertex> ,peano::stacks::CellArrayStack<exahype::Cell> ,exahype::adapters::Correction>;
    template class Grid<exahype::Vertex,exahype::Cell,exahype::State, peano::stacks::VertexSTDStack<  exahype::Vertex> ,peano::stacks::CellSTDStack<  exahype::Cell> ,exahype::adapters::Correction>;
  }
}

#include "peano/grid/Grid.cpph"
