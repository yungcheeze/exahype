#include "ExplicitEulerForHeatEquation/repositories/Repository.h"
#include "ExplicitEulerForHeatEquation/records/RepositoryState.h"

#include "ExplicitEulerForHeatEquation/State.h"
#include "ExplicitEulerForHeatEquation/Vertex.h"
#include "ExplicitEulerForHeatEquation/Cell.h"

#include "peano/grid/Grid.h"

#include "peano/stacks/CellArrayStack.h"
#include "peano/stacks/CellSTDStack.h"

#include "peano/stacks/VertexArrayStack.h"
#include "peano/stacks/VertexSTDStack.h"

 #include "ExplicitEulerForHeatEquation/adapters/CreateGrid.h" 
 #include "ExplicitEulerForHeatEquation/adapters/TimeStep.h" 
 #include "ExplicitEulerForHeatEquation/adapters/CreateGridAndPlot.h" 
 #include "ExplicitEulerForHeatEquation/adapters/TimeStepAndPlot.h" 


namespace peano {
  namespace grid {
    template class Grid<myproject::Vertex,myproject::Cell,myproject::State, peano::stacks::VertexArrayStack<myproject::Vertex> ,peano::stacks::CellArrayStack<myproject::Cell> ,myproject::adapters::CreateGrid>;
    template class Grid<myproject::Vertex,myproject::Cell,myproject::State, peano::stacks::VertexSTDStack<  myproject::Vertex> ,peano::stacks::CellSTDStack<  myproject::Cell> ,myproject::adapters::CreateGrid>;
  }
}

#include "peano/grid/Grid.cpph"
