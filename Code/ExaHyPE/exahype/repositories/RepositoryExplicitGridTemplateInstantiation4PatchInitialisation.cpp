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

 #include "exahype/adapters/InitialGrid.h" 
 #include "exahype/adapters/GridExport.h" 
 #include "exahype/adapters/PatchInitialisation.h" 
 #include "exahype/adapters/PatchInitialisationAndExport.h" 
 #include "exahype/adapters/FaceDataExchange.h" 
 #include "exahype/adapters/InitialConditionAndGlobalTimeStepComputation.h" 
 #include "exahype/adapters/InitialConditionAndExportAndGlobalTimeStepComputation.h" 
 #include "exahype/adapters/PredictorAndGlobalTimeStepComputation.h" 
 #include "exahype/adapters/CorrectorAndPredictorAndGlobalTimeStepComputation.h" 
 #include "exahype/adapters/CorrectorAndPredictorAndGlobalTimeStepComputationAndExport.h" 


namespace peano {
  namespace grid {
    template class Grid<exahype::Vertex,exahype::Cell,exahype::State, peano::stacks::VertexArrayStack<exahype::Vertex> ,peano::stacks::CellArrayStack<exahype::Cell> ,exahype::adapters::PatchInitialisation>;
    template class Grid<exahype::Vertex,exahype::Cell,exahype::State, peano::stacks::VertexSTDStack<  exahype::Vertex> ,peano::stacks::CellSTDStack<  exahype::Cell> ,exahype::adapters::PatchInitialisation>;
  }
}

#include "peano/grid/Grid.cpph"
