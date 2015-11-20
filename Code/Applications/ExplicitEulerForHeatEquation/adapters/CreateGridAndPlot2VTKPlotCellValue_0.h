// This file is part of the Peano project. For conditions of distribution and 
// use, please see the copyright notice at www.peano-framework.org
#ifndef MYPROJECT_ADAPTERS_CreateGridAndPlot2VTKPlotCellValue_0_H_
#define MYPROJECT_ADAPTERS_CreateGridAndPlot2VTKPlotCellValue_0_H_


#include "tarch/la/Vector.h"
#include "tarch/la/VectorCompare.h"
#include "tarch/logging/Log.h"
#include "tarch/multicore/MulticoreDefinitions.h"
#include "tarch/plotter/griddata/unstructured/vtk/VTKTextFileWriter.h"

#include "peano/MappingSpecification.h"
#include "peano/CommunicationSpecification.h"
#include "peano/grid/VertexEnumerator.h"

#include "ExplicitEulerForHeatEquation/Vertex.h"
#include "ExplicitEulerForHeatEquation/Cell.h"
#include "ExplicitEulerForHeatEquation/State.h"

#include <map>


namespace myproject {
      namespace adapters {
        class CreateGridAndPlot2VTKPlotCellValue_0;
      } 
}


/**
 * This is an adapter plotting a vtk grid file. Please set
 *
 * epsilon   filename
 * getEpsilon   getter for cell value (without brackets)
 * eps   Legend title
 *
 * @author Tobias Weinzierl
 * @version $Revision: 1.10 $
 */
class myproject::adapters::CreateGridAndPlot2VTKPlotCellValue_0 {
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
    
    tarch::plotter::griddata::unstructured::vtk::VTKTextFileWriter*                 _vtkWriter;
    tarch::plotter::griddata::unstructured::UnstructuredGridWriter::VertexWriter*   _vertexWriter;
    tarch::plotter::griddata::unstructured::UnstructuredGridWriter::CellWriter*     _cellWriter;
    
    tarch::plotter::griddata::Writer::CellDataWriter*                               _cellValueWriter;
    
    static int _snapshotCounter;
    
    void plotVertex(
      const myproject::Vertex&                 fineGridVertex,
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

    CreateGridAndPlot2VTKPlotCellValue_0();

    #if defined(SharedMemoryParallelisation)
    CreateGridAndPlot2VTKPlotCellValue_0(const CreateGridAndPlot2VTKPlotCellValue_0& masterThread);
    #endif

    virtual ~CreateGridAndPlot2VTKPlotCellValue_0();
  
    #if defined(SharedMemoryParallelisation)
    void mergeWithWorkerThread(const CreateGridAndPlot2VTKPlotCellValue_0& workerThread);
    #endif

    void createInnerVertex(
      myproject::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      myproject::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      myproject::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
    );


    void createBoundaryVertex(
      myproject::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      myproject::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      myproject::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
    );


    void createHangingVertex(
      myproject::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      myproject::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      myproject::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
    );


    void destroyHangingVertex(
      const myproject::Vertex&   fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      myproject::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      myproject::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
    );


    void destroyVertex(
      const myproject::Vertex&   fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      myproject::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      myproject::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
    );


    void createCell(
      myproject::Cell&                 fineGridCell,
      myproject::Vertex * const         fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      myproject::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      myproject::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
    );


    void destroyCell(
      const myproject::Cell&           fineGridCell,
      myproject::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      myproject::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      myproject::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
    );
        
    #ifdef Parallel
    void mergeWithNeighbour(
      myproject::Vertex&  vertex,
      const myproject::Vertex&  neighbour,
      int                                           fromRank,
      const tarch::la::Vector<DIMENSIONS,double>&   fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&   fineGridH,
      int                                           level
    );

    void prepareSendToNeighbour(
      myproject::Vertex&  vertex,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
    );

    void prepareCopyToRemoteNode(
      myproject::Vertex&  localVertex,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
    );

    void prepareCopyToRemoteNode(
      myproject::Cell&  localCell,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   cellCentre,
      const tarch::la::Vector<DIMENSIONS,double>&   cellSize,
      int                                           level
    );

    void mergeWithRemoteDataDueToForkOrJoin(
      myproject::Vertex&  localVertex,
      const myproject::Vertex&  masterOrWorkerVertex,
      int                                       fromRank,
      const tarch::la::Vector<DIMENSIONS,double>&  x,
      const tarch::la::Vector<DIMENSIONS,double>&  h,
      int                                       level
    );

    void mergeWithRemoteDataDueToForkOrJoin(
      myproject::Cell&  localCell,
      const myproject::Cell&  masterOrWorkerCell,
      int                                       fromRank,
      const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
      const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
      int                                       level
    );

    bool prepareSendToWorker(
      myproject::Cell&                 fineGridCell,
      myproject::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      myproject::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      myproject::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell,
      int                                                                  worker
    );

    void prepareSendToMaster(
      myproject::Cell&                       localCell,
      myproject::Vertex *                    vertices,
      const peano::grid::VertexEnumerator&       verticesEnumerator, 
      const myproject::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&       coarseGridVerticesEnumerator,
      const myproject::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&   fineGridPositionOfCell
    );

    void mergeWithMaster(
      const myproject::Cell&           workerGridCell,
      myproject::Vertex * const        workerGridVertices,
      const peano::grid::VertexEnumerator& workerEnumerator,
      myproject::Cell&                 fineGridCell,
      myproject::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      myproject::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      myproject::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell,
      int                                                                  worker,
      const myproject::State&          workerState,
      myproject::State&                masterState
    );


    void receiveDataFromMaster(
      myproject::Cell&                        receivedCell, 
      myproject::Vertex *                     receivedVertices,
      const peano::grid::VertexEnumerator&        receivedVerticesEnumerator,
      myproject::Vertex * const               receivedCoarseGridVertices,
      const peano::grid::VertexEnumerator&        receivedCoarseGridVerticesEnumerator,
      myproject::Cell&                        receivedCoarseGridCell,
      myproject::Vertex * const               workersCoarseGridVertices,
      const peano::grid::VertexEnumerator&        workersCoarseGridVerticesEnumerator,
      myproject::Cell&                        workersCoarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&    fineGridPositionOfCell
    );


    void mergeWithWorker(
      myproject::Cell&           localCell, 
      const myproject::Cell&     receivedMasterCell,
      const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
      const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
      int                                          level
    );


    void mergeWithWorker(
      myproject::Vertex&        localVertex,
      const myproject::Vertex&  receivedMasterVertex,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
    );
    #endif


    void touchVertexFirstTime(
      myproject::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      myproject::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      myproject::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
    );


    void touchVertexLastTime(
      myproject::Vertex&         fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      myproject::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      myproject::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
    );
    

    void enterCell(
      myproject::Cell&                 fineGridCell,
      myproject::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      myproject::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      myproject::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
    );


    void leaveCell(
      myproject::Cell&                          fineGridCell,
      myproject::Vertex * const                 fineGridVertices,
      const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
      myproject::Vertex * const                 coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      myproject::Cell&                          coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&      fineGridPositionOfCell
    );


    void beginIteration(
      myproject::State&  solverState
    );


    void endIteration(
      myproject::State&  solverState
    );

    void descend(
      myproject::Cell * const          fineGridCells,
      myproject::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      myproject::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      myproject::Cell&                 coarseGridCell
    );


    void ascend(
      myproject::Cell * const    fineGridCells,
      myproject::Vertex * const  fineGridVertices,
      const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
      myproject::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      myproject::Cell&           coarseGridCell
    );    
};


#endif
