#include "AdvectionDG0/repositories/Repository.h"
#include "AdvectionDG0/records/RepositoryState.h"

#include "AdvectionDG0/State.h"
#include "AdvectionDG0/Vertex.h"
#include "AdvectionDG0/Cell.h"

#include "peano/grid/Grid.h"

#include "peano/stacks/CellArrayStack.h"
#include "peano/stacks/CellSTDStack.h"

#include "peano/stacks/VertexArrayStack.h"
#include "peano/stacks/VertexSTDStack.h"

 #include "AdvectionDG0/adapters/CreateGrid.h" 
 #include "AdvectionDG0/adapters/TimeStep.h" 


namespace peano {
  namespace grid {
    template class Grid<exahype::Vertex,exahype::Cell,exahype::State, peano::stacks::VertexArrayStack<exahype::Vertex> ,peano::stacks::CellArrayStack<exahype::Cell> ,exahype::adapters::TimeStep>;
    template class Grid<exahype::Vertex,exahype::Cell,exahype::State, peano::stacks::VertexSTDStack<  exahype::Vertex> ,peano::stacks::CellSTDStack<  exahype::Cell> ,exahype::adapters::TimeStep>;
  }
}

#include "peano/grid/Grid.cpph"
