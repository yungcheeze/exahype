#include "ExaHyPeDataStructure/repositories/Repository.h"
#include "ExaHyPeDataStructure/records/RepositoryState.h"

#include "ExaHyPeDataStructure/State.h"
#include "ExaHyPeDataStructure/Vertex.h"
#include "ExaHyPeDataStructure/Cell.h"

#include "peano/grid/Grid.h"

#include "peano/stacks/CellArrayStack.h"
#include "peano/stacks/CellSTDStack.h"

#include "peano/stacks/VertexArrayStack.h"
#include "peano/stacks/VertexSTDStack.h"

 #include "ExaHyPeDataStructure/adapters/CreateGrid.h" 
 #include "ExaHyPeDataStructure/adapters/TimeStep.h" 


namespace peano {
  namespace grid {
    template class Grid<exahype::Vertex,exahype::Cell,exahype::State, peano::stacks::VertexArrayStack<exahype::Vertex> ,peano::stacks::CellArrayStack<exahype::Cell> ,exahype::adapters::TimeStep>;
    template class Grid<exahype::Vertex,exahype::Cell,exahype::State, peano::stacks::VertexSTDStack<  exahype::Vertex> ,peano::stacks::CellSTDStack<  exahype::Cell> ,exahype::adapters::TimeStep>;
  }
}

#include "peano/grid/Grid.cpph"
