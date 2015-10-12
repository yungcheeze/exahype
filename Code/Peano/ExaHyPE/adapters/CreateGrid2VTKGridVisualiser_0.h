// This file is part of the Peano project. For conditions of distribution and 
// use, please see the copyright notice at www.peano-framework.org
#ifndef EXAHYPE_ADAPTERS_CreateGrid2VTKGridVisualiser_0_H_
#define EXAHYPE_ADAPTERS_CreateGrid2VTKGridVisualiser_0_H_


#include "tarch/la/Vector.h"
#include "tarch/la/VectorCompare.h"
#include "tarch/logging/Log.h"
#include "tarch/multicore/MulticoreDefinitions.h"
#include "tarch/plotter/griddata/unstructured/vtk/VTKBinaryFileWriter.h"

#include "peano/MappingSpecification.h"
#include "peano/CommunicationSpecification.h"
#include "peano/grid/VertexEnumerator.h"

#include "ExaHyPE/Vertex.h"
#include "ExaHyPE/Cell.h"
#include "ExaHyPE/State.h"

#include <map>


namespace ExaHyPE {
      namespace adapters {
        class CreateGrid2VTKGridVisualiser_0;
      } 
}


/**
 * This is an adapter plotting a vtk grid file. Please set
 *
 * grid   filename
 *
 * @author Tobias Weinzierl
 * @version $Revision: 1.10 $
 */
class ExaHyPE::adapters::CreateGrid2VTKGridVisualiser_0 {
  private:
    /**
     * One big map mapping vertices to indices. The procedure using this map is 
     * straightforward. Whenever we encounter a vertex, the object does a 
     * lookup whether this vertex already has been plotted. If not, it plots it 
     * and adds an entry.
     * 
     * @see plotVertex(const tarch::la::Vector<DIMENSIONS,double>&  x)
     */
    static std::map<tarch::la::Vector<DIMENSIONS,double> , int, tarch::la::VectorCompare<DIMENSIONS> >  _vertex2IndexMap;
    
    tarch::plotter::griddata::unstructured::vtk::VTKBinaryFileWriter*               _vtkWriter;
    tarch::plotter::griddata::unstructured::UnstructuredGridWriter::VertexWriter*   _vertexWriter;
    tarch::plotter::griddata::unstructured::UnstructuredGridWriter::CellWriter*     _cellWriter;
    
    tarch::plotter::griddata::Writer::VertexDataWriter*                             _vertexTypeWriter;
    tarch::plotter::griddata::Writer::VertexDataWriter*                             _vertexRefinementControlWriter;
    tarch::plotter::griddata::Writer::VertexDataWriter*                             _vertexAdjacentCellsHeight;

    tarch::plotter::griddata::Writer::CellDataWriter*                               _cellStateWriter;
    
    static int _snapshotCounter;
    
    void plotVertex(
      const ExaHyPE::Vertex&                 fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&  fineGridX
    );
  public:
    static peano::MappingSpecification   touchVertexLastTimeSpecification();
    static peano::MappingSpecification   touchVertexFirstTimeSpecification();
    static peano::MappingSpecification   enterCellSpecification();
    static peano::MappingSpecification   leaveCellSpecification();
    static peano::MappingSpecification   ascendSpecification();
    static peano::MappingSpecification   descendSpecification();

    static peano::CommunicationSpecification   communicationSpecification();

    CreateGrid2VTKGridVisualiser_0();

    #if defined(SharedMemoryParallelisation)
    CreateGrid2VTKGridVisualiser_0(const CreateGrid2VTKGridVisualiser_0& masterThread);
    #endif

    virtual ~CreateGrid2VTKGridVisualiser_0();
  
    #if defined(SharedMemoryParallelisation)
    void mergeWithWorkerThread(const CreateGrid2VTKGridVisualiser_0& workerThread);
    #endif

    void createInnerVertex(
      ExaHyPE::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      ExaHyPE::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      ExaHyPE::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
    );


    void createBoundaryVertex(
      ExaHyPE::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      ExaHyPE::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      ExaHyPE::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
    );


    void createHangingVertex(
      ExaHyPE::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      ExaHyPE::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      ExaHyPE::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
    );


    void destroyHangingVertex(
      const ExaHyPE::Vertex&   fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      ExaHyPE::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      ExaHyPE::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
    );


    void destroyVertex(
      const ExaHyPE::Vertex&   fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      ExaHyPE::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      ExaHyPE::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
    );


    void createCell(
      ExaHyPE::Cell&                 fineGridCell,
      ExaHyPE::Vertex * const         fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      ExaHyPE::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      ExaHyPE::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
    );


    void destroyCell(
      const ExaHyPE::Cell&           fineGridCell,
      ExaHyPE::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      ExaHyPE::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      ExaHyPE::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
    );
        
    #ifdef Parallel
    void mergeWithNeighbour(
      ExaHyPE::Vertex&  vertex,
      const ExaHyPE::Vertex&  neighbour,
      int                                           fromRank,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
    );

    void prepareSendToNeighbour(
      ExaHyPE::Vertex&  vertex,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
    );

    void prepareCopyToRemoteNode(
      ExaHyPE::Vertex&  localVertex,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
    );

    void prepareCopyToRemoteNode(
      ExaHyPE::Cell&  localCell,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   cellCentre,
      const tarch::la::Vector<DIMENSIONS,double>&   cellSize,
      int                                           level
    );

    void mergeWithRemoteDataDueToForkOrJoin(
      ExaHyPE::Vertex&  localVertex,
      const ExaHyPE::Vertex&  masterOrWorkerVertex,
      int                                       fromRank,
      const tarch::la::Vector<DIMENSIONS,double>&  x,
      const tarch::la::Vector<DIMENSIONS,double>&  h,
      int                                       level
    );

    void mergeWithRemoteDataDueToForkOrJoin(
      ExaHyPE::Cell&  localCell,
      const ExaHyPE::Cell&  masterOrWorkerCell,
      int                                       fromRank,
      const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
      const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
      int                                       level
    );

    bool prepareSendToWorker(
      ExaHyPE::Cell&                 fineGridCell,
      ExaHyPE::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      ExaHyPE::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      ExaHyPE::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell,
      int                                                                  worker
    );

    void prepareSendToMaster(
      ExaHyPE::Cell&                       localCell,
      ExaHyPE::Vertex *                    vertices,
      const peano::grid::VertexEnumerator&       verticesEnumerator, 
      const ExaHyPE::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&       coarseGridVerticesEnumerator,
      const ExaHyPE::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&   fineGridPositionOfCell
    );

    void mergeWithMaster(
      const ExaHyPE::Cell&           workerGridCell,
      ExaHyPE::Vertex * const        workerGridVertices,
      const peano::grid::VertexEnumerator& workerEnumerator,
      ExaHyPE::Cell&                 fineGridCell,
      ExaHyPE::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      ExaHyPE::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      ExaHyPE::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell,
      int                                                                  worker,
      const ExaHyPE::State&          workerState,
      ExaHyPE::State&                masterState
    );


    void receiveDataFromMaster(
      ExaHyPE::Cell&                        receivedCell, 
      ExaHyPE::Vertex *                     receivedVertices,
      const peano::grid::VertexEnumerator&        receivedVerticesEnumerator,
      ExaHyPE::Vertex * const               receivedCoarseGridVertices,
      const peano::grid::VertexEnumerator&        receivedCoarseGridVerticesEnumerator,
      ExaHyPE::Cell&                        receivedCoarseGridCell,
      ExaHyPE::Vertex * const               workersCoarseGridVertices,
      const peano::grid::VertexEnumerator&        workersCoarseGridVerticesEnumerator,
      ExaHyPE::Cell&                        workersCoarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&    fineGridPositionOfCell
    );


    void mergeWithWorker(
      ExaHyPE::Cell&           localCell, 
      const ExaHyPE::Cell&     receivedMasterCell,
      const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
      const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
      int                                          level
    );


    void mergeWithWorker(
      ExaHyPE::Vertex&        localVertex,
      const ExaHyPE::Vertex&  receivedMasterVertex,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
    );
    #endif


    void touchVertexFirstTime(
      ExaHyPE::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      ExaHyPE::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      ExaHyPE::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
    );


    void touchVertexLastTime(
      ExaHyPE::Vertex&         fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      ExaHyPE::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      ExaHyPE::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
    );
    

    void enterCell(
      ExaHyPE::Cell&                 fineGridCell,
      ExaHyPE::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      ExaHyPE::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      ExaHyPE::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
    );


    void leaveCell(
      ExaHyPE::Cell&                          fineGridCell,
      ExaHyPE::Vertex * const                 fineGridVertices,
      const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
      ExaHyPE::Vertex * const                 coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      ExaHyPE::Cell&                          coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&      fineGridPositionOfCell
    );


    void beginIteration(
      ExaHyPE::State&  solverState
    );


    void endIteration(
      ExaHyPE::State&  solverState
    );

    void descend(
      ExaHyPE::Cell * const          fineGridCells,
      ExaHyPE::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      ExaHyPE::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      ExaHyPE::Cell&                 coarseGridCell
    );


    void ascend(
      ExaHyPE::Cell * const    fineGridCells,
      ExaHyPE::Vertex * const  fineGridVertices,
      const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
      ExaHyPE::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      ExaHyPE::Cell&           coarseGridCell
    );    
};


#endif
