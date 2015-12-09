#include "EulerFlow3d/repositories/Repository.h"
#include "EulerFlow3d/records/RepositoryState.h"

#include "EulerFlow3d/State.h"
#include "EulerFlow3d/Vertex.h"
#include "EulerFlow3d/Cell.h"

#include "peano/grid/Grid.h"

#include "peano/stacks/CellArrayStack.h"
#include "peano/stacks/CellSTDStack.h"

#include "peano/stacks/VertexArrayStack.h"
#include "peano/stacks/VertexSTDStack.h"

 #include "EulerFlow3d/adapters/InitialGrid.h" 
 #include "EulerFlow3d/adapters/GridExport.h" 
 #include "EulerFlow3d/adapters/PatchInitialisation.h" 
 #include "EulerFlow3d/adapters/PatchInitialisationAndExport.h" 
 #include "EulerFlow3d/adapters/InitialCondition.h" 
 #include "EulerFlow3d/adapters/InitialConditionAndExport.h" 
 #include "EulerFlow3d/adapters/GlobalTimeStepComputation.h" 
 #include "EulerFlow3d/adapters/Predictor.h" 
 #include "EulerFlow3d/adapters/Corrector.h" 
 #include "EulerFlow3d/adapters/CorrectorAndExport.h" 
 #include "EulerFlow3d/adapters/SolutionExport.h" 


namespace peano {
  namespace grid {
    template class Grid<exahype::Vertex,exahype::Cell,exahype::State, peano::stacks::VertexArrayStack<exahype::Vertex> ,peano::stacks::CellArrayStack<exahype::Cell> ,exahype::adapters::SolutionExport>;
    template class Grid<exahype::Vertex,exahype::Cell,exahype::State, peano::stacks::VertexSTDStack<  exahype::Vertex> ,peano::stacks::CellSTDStack<  exahype::Cell> ,exahype::adapters::SolutionExport>;
  }
}

#include "peano/grid/Grid.cpph"
