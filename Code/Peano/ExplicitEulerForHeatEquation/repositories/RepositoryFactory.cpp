#include "ExplicitEulerForHeatEquation/repositories/RepositoryFactory.h"

#include "ExplicitEulerForHeatEquation/repositories/RepositoryArrayStack.h"
#include "ExplicitEulerForHeatEquation/repositories/RepositorySTDStack.h"

#include "ExplicitEulerForHeatEquation/records/RepositoryState.h"

#ifdef Parallel
#include "tarch/parallel/NodePool.h"
#include "peano/parallel/Partitioner.h"
#endif


myproject::repositories::RepositoryFactory::RepositoryFactory() {
  #ifdef Parallel
  peano::parallel::Partitioner::initDatatypes();

  myproject::State::initDatatype();
  myproject::Vertex::initDatatype();
  myproject::Cell::initDatatype();

  if (myproject::records::RepositoryState::Datatype==0) {
    myproject::records::RepositoryState::initDatatype();
  }
  #endif
}


myproject::repositories::RepositoryFactory::~RepositoryFactory() {
}


void myproject::repositories::RepositoryFactory::shutdownAllParallelDatatypes() {
  #ifdef Parallel
  peano::parallel::Partitioner::shutdownDatatypes();

  myproject::State::shutdownDatatype();
  myproject::Vertex::shutdownDatatype();
  myproject::Cell::shutdownDatatype();

  if (myproject::records::RepositoryState::Datatype!=0) {
    myproject::records::RepositoryState::shutdownDatatype();
    myproject::records::RepositoryState::Datatype = 0;
  }
  #endif
}


myproject::repositories::RepositoryFactory& myproject::repositories::RepositoryFactory::getInstance() {
  static myproject::repositories::RepositoryFactory singleton;
  return singleton;
}

    
myproject::repositories::Repository* 
myproject::repositories::RepositoryFactory::createWithArrayStackImplementation(
  peano::geometry::Geometry&                   geometry,
  const tarch::la::Vector<DIMENSIONS,double>&  domainSize,
  const tarch::la::Vector<DIMENSIONS,double>&  computationalDomainOffset,
  int                                          maxCellStackSize,    
  int                                          maxVertexStackSize,    
  int                                          maxTemporaryVertexStackSize    
) {
  #ifdef Parallel
  if (!tarch::parallel::Node::getInstance().isGlobalMaster()) {
    return new myproject::repositories::RepositoryArrayStack(geometry, domainSize, computationalDomainOffset,maxCellStackSize,maxVertexStackSize,maxTemporaryVertexStackSize);
  }
  else
  #endif
  return new myproject::repositories::RepositoryArrayStack(geometry, domainSize, computationalDomainOffset,maxCellStackSize,maxVertexStackSize,maxTemporaryVertexStackSize);
}    


myproject::repositories::Repository* 
myproject::repositories::RepositoryFactory::createWithSTDStackImplementation(
  peano::geometry::Geometry&                   geometry,
  const tarch::la::Vector<DIMENSIONS,double>&  domainSize,
  const tarch::la::Vector<DIMENSIONS,double>&  computationalDomainOffset
) {
  #ifdef Parallel
  if (!tarch::parallel::Node::getInstance().isGlobalMaster()) {
    return new myproject::repositories::RepositorySTDStack(geometry);
  }
  else
  #endif
  return new myproject::repositories::RepositorySTDStack(geometry, domainSize, computationalDomainOffset);
}
