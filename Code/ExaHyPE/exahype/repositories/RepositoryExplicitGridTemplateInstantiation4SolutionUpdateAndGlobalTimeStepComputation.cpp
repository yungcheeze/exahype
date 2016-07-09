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

 #include "exahype/adapters/RegularMesh.h" 
 #include "exahype/adapters/AugmentedAMRGrid.h" 
 #include "exahype/adapters/PlotAugmentedAMRGrid.h" 
 #include "exahype/adapters/SolutionUpdateAndGlobalTimeStepComputation.h" 
 #include "exahype/adapters/PredictorAndPlotAndGlobalTimeStepComputation.h" 
 #include "exahype/adapters/PredictorAndGlobalTimeStepComputation.h" 
 #include "exahype/adapters/ADERDGTimeStep.h" 
 #include "exahype/adapters/ADERDGTimeStepAndPlot.h" 
 #include "exahype/adapters/GlobalTimeStepComputation.h" 
 #include "exahype/adapters/PlotAndGlobalTimeStepComputation.h" 
 #include "exahype/adapters/FaceDataExchange.h" 
 #include "exahype/adapters/Predictor.h" 
 #include "exahype/adapters/PredictorRerun.h" 
 #include "exahype/adapters/Corrector.h" 


namespace peano {
  namespace grid {
    template class Grid<exahype::Vertex,exahype::Cell,exahype::State, peano::stacks::VertexArrayStack<exahype::Vertex> ,peano::stacks::CellArrayStack<exahype::Cell> ,exahype::adapters::SolutionUpdateAndGlobalTimeStepComputation>;
    template class Grid<exahype::Vertex,exahype::Cell,exahype::State, peano::stacks::VertexSTDStack<  exahype::Vertex> ,peano::stacks::CellSTDStack<  exahype::Cell> ,exahype::adapters::SolutionUpdateAndGlobalTimeStepComputation>;
  }
}

#include "peano/grid/Grid.cpph"
