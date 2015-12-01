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

 #include "EulerFlow3d/adapters/CreateGrid.h" 
 #include "EulerFlow3d/adapters/PlotGrid.h" 
 #include "EulerFlow3d/adapters/InitCells.h" 
 #include "EulerFlow3d/adapters/InitCellData.h" 
 #include "EulerFlow3d/adapters/TimeStep.h" 
 #include "EulerFlow3d/adapters/TimeStepAndPlot.h" 


namespace peano {
  namespace grid {
    template class Grid<exahype::Vertex,exahype::Cell,exahype::State, peano::stacks::VertexArrayStack<exahype::Vertex> ,peano::stacks::CellArrayStack<exahype::Cell> ,exahype::adapters::InitCells>;
    template class Grid<exahype::Vertex,exahype::Cell,exahype::State, peano::stacks::VertexSTDStack<  exahype::Vertex> ,peano::stacks::CellSTDStack<  exahype::Cell> ,exahype::adapters::InitCells>;
  }
}

#include "peano/grid/Grid.cpph"
