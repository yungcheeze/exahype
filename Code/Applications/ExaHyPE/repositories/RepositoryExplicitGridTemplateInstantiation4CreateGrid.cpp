#include "ExaHyPE/repositories/Repository.h"
#include "ExaHyPE/records/RepositoryState.h"

#include "ExaHyPE/State.h"
#include "ExaHyPE/Vertex.h"
#include "ExaHyPE/Cell.h"

#include "peano/grid/Grid.h"

#include "peano/stacks/CellArrayStack.h"
#include "peano/stacks/CellSTDStack.h"

#include "peano/stacks/VertexArrayStack.h"
#include "peano/stacks/VertexSTDStack.h"

 #include "ExaHyPE/adapters/CreateGrid.h" 


namespace peano {
  namespace grid {
    template class Grid<ExaHyPE::Vertex,ExaHyPE::Cell,ExaHyPE::State, peano::stacks::VertexArrayStack<ExaHyPE::Vertex> ,peano::stacks::CellArrayStack<ExaHyPE::Cell> ,ExaHyPE::adapters::CreateGrid>;
    template class Grid<ExaHyPE::Vertex,ExaHyPE::Cell,ExaHyPE::State, peano::stacks::VertexSTDStack<  ExaHyPE::Vertex> ,peano::stacks::CellSTDStack<  ExaHyPE::Cell> ,ExaHyPE::adapters::CreateGrid>;
  }
}

#include "peano/grid/Grid.cpph"
