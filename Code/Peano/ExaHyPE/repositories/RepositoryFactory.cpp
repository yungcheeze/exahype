#include "ExaHyPE/repositories/RepositoryFactory.h"

#include "ExaHyPE/repositories/RepositoryArrayStack.h"
#include "ExaHyPE/repositories/RepositorySTDStack.h"

#include "ExaHyPE/records/RepositoryState.h"

#ifdef Parallel
#include "tarch/parallel/NodePool.h"
#include "peano/parallel/Partitioner.h"
#endif


ExaHyPE::repositories::RepositoryFactory::RepositoryFactory() {
  #ifdef Parallel
  peano::parallel::Partitioner::initDatatypes();

  ExaHyPE::State::initDatatype();
  ExaHyPE::Vertex::initDatatype();
  ExaHyPE::Cell::initDatatype();

  if (ExaHyPE::records::RepositoryState::Datatype==0) {
    ExaHyPE::records::RepositoryState::initDatatype();
  }
  #endif
}


ExaHyPE::repositories::RepositoryFactory::~RepositoryFactory() {
}


void ExaHyPE::repositories::RepositoryFactory::shutdownAllParallelDatatypes() {
  #ifdef Parallel
  peano::parallel::Partitioner::shutdownDatatypes();

  ExaHyPE::State::shutdownDatatype();
  ExaHyPE::Vertex::shutdownDatatype();
  ExaHyPE::Cell::shutdownDatatype();

  if (ExaHyPE::records::RepositoryState::Datatype!=0) {
    ExaHyPE::records::RepositoryState::shutdownDatatype();
    ExaHyPE::records::RepositoryState::Datatype = 0;
  }
  #endif
}


ExaHyPE::repositories::RepositoryFactory& ExaHyPE::repositories::RepositoryFactory::getInstance() {
  static ExaHyPE::repositories::RepositoryFactory singleton;
  return singleton;
}

    
ExaHyPE::repositories::Repository* 
ExaHyPE::repositories::RepositoryFactory::createWithArrayStackImplementation(
  peano::geometry::Geometry&                   geometry,
  const tarch::la::Vector<DIMENSIONS,double>&  domainSize,
  const tarch::la::Vector<DIMENSIONS,double>&  computationalDomainOffset,
  int                                          maxCellStackSize,    
  int                                          maxVertexStackSize,    
  int                                          maxTemporaryVertexStackSize    
) {
  #ifdef Parallel
  if (!tarch::parallel::Node::getInstance().isGlobalMaster()) {
    return new ExaHyPE::repositories::RepositoryArrayStack(geometry, domainSize, computationalDomainOffset,maxCellStackSize,maxVertexStackSize,maxTemporaryVertexStackSize);
  }
  else
  #endif
  return new ExaHyPE::repositories::RepositoryArrayStack(geometry, domainSize, computationalDomainOffset,maxCellStackSize,maxVertexStackSize,maxTemporaryVertexStackSize);
}    


ExaHyPE::repositories::Repository* 
ExaHyPE::repositories::RepositoryFactory::createWithSTDStackImplementation(
  peano::geometry::Geometry&                   geometry,
  const tarch::la::Vector<DIMENSIONS,double>&  domainSize,
  const tarch::la::Vector<DIMENSIONS,double>&  computationalDomainOffset
) {
  #ifdef Parallel
  if (!tarch::parallel::Node::getInstance().isGlobalMaster()) {
    return new ExaHyPE::repositories::RepositorySTDStack(geometry);
  }
  else
  #endif
  return new ExaHyPE::repositories::RepositorySTDStack(geometry, domainSize, computationalDomainOffset);
}
