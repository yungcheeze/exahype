#include "EulerFlow/repositories/Repository.h"
#include "EulerFlow/records/RepositoryState.h"

#include "EulerFlow/State.h"
#include "EulerFlow/Vertex.h"
#include "EulerFlow/Cell.h"

#include "peano/grid/Grid.h"

#include "peano/stacks/CellArrayStack.h"
#include "peano/stacks/CellSTDStack.h"

#include "peano/stacks/VertexArrayStack.h"
#include "peano/stacks/VertexSTDStack.h"

 #include "EulerFlow/adapters/InitialGrid.h" 
 #include "EulerFlow/adapters/GridExport.h" 
 #include "EulerFlow/adapters/PatchInitialisation.h" 
 #include "EulerFlow/adapters/PatchInitialisationAndExport.h" 
 #include "EulerFlow/adapters/FaceDataExchange.h" 
 #include "EulerFlow/adapters/InitialConditionAndGlobalTimeStepComputation.h" 
 #include "EulerFlow/adapters/InitialConditionAndExportAndGlobalTimeStepComputation.h" 
 #include "EulerFlow/adapters/PredictorAndGlobalTimeStepComputation.h" 
 #include "EulerFlow/adapters/CorrectorAndPredictorAndGlobalTimeStepComputation.h" 
 #include "EulerFlow/adapters/CorrectorAndPredictorAndGlobalTimeStepComputationAndExport.h" 


namespace peano {
  namespace grid {
    template class Grid<exahype::Vertex,exahype::Cell,exahype::State, peano::stacks::VertexArrayStack<exahype::Vertex> ,peano::stacks::CellArrayStack<exahype::Cell> ,exahype::adapters::InitialGrid>;
    template class Grid<exahype::Vertex,exahype::Cell,exahype::State, peano::stacks::VertexSTDStack<  exahype::Vertex> ,peano::stacks::CellSTDStack<  exahype::Cell> ,exahype::adapters::InitialGrid>;
  }
}

#include "peano/grid/Grid.cpph"
