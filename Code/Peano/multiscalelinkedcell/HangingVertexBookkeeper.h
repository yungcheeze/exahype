// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www.peano-framework.org
#ifndef _MULISCALE_LINKED_CELL_HANGING_VERTEX_BOOKKEEPER_H_
#define _MULISCALE_LINKED_CELL_HANGING_VERTEX_BOOKKEEPER_H_


#include "tarch/logging/Log.h"

#include "tarch/la/Vector.h"
#include "tarch/la/VectorCompare.h"

#include "peano/utils/Globals.h"
#include "peano/grid/VertexEnumerator.h"

#include <map>

#define CompilerDoesNotSupportC11ContainerErase

namespace multiscalelinkedcell {
  class HangingVertexBookkeeper;

  /**
   * Return the indices of the cells surrounding a cell
   *
   * Throw in the 2^d * 2^d vector returned by an autogenerated readXXX
   * operation on the vertices and extract the 3^d-1 cells surrounding the
   * central cell (the one additional entry contains the index of the cell.
   * The result is enumerated lexicographically:
   *
   * @image html getIndicesAroundCell.png
   *
   * The image above illustrates the behaviour for d=2. You hand in cell
   * indices belonging to the coloured cell. The operation returns a
   * vector of size nine. The very first entry is the left bottom neighbour,
   * the second is the bottom neighbour, and so forth. The fifth entry (index
   * four) is the cell's index itself, and so forth.
   *
   * In 3d, the enumeration follows straightforwardly: We first run along
   * the x axis, then along the y axis, then along the z axis. The fifth
   * entry then is the cell in front of the current cell, the eleventh entry
   * is the bottom neigbhour, and the 14th entry is the cell index itself.
   */
  tarch::la::Vector<THREE_POWER_D,int> getIndicesAroundCell(
    const tarch::la::Vector<TWO_POWER_D_TIMES_TWO_POWER_D,int>&  indices
  );
}


/**
 * Singleon keeping track of the adjacency data of the hanging nodes.
 *
 * To get the bookkeeper's instance, i.e. to work with it, you call
 * \code
 multiscalelinkedcell::HangingVertexBookkeeper::getInstance().
   \endcode
 *
 * @author Kristof Unterweger, Tobias Weinzierl
 */
class multiscalelinkedcell::HangingVertexBookkeeper {
  private:
    static tarch::logging::Log  _log;

    struct HangingVertexIdentifier {
      tarch::la::Vector<TWO_POWER_D,int>  indicesOfAdjacentCells;
      bool                                usedInLastTraversal;
    };

    typedef std::map< tarch::la::Vector<DIMENSIONS+1,double >, HangingVertexIdentifier, tarch::la::VectorCompare<DIMENSIONS+1> >  VertexMap;

    VertexMap _vertexMap;

    HangingVertexBookkeeper();

    tarch::la::Vector<DIMENSIONS+1,double > getKey(
      const tarch::la::Vector<DIMENSIONS,double>&  x,
      int                                          level
    ) const;

  public:
    /**
     * Every index greater is valid
     */
    static const int InvalidAdjacencyIndex;
    static const int RemoteAdjacencyIndex;
    static const int DomainBoundaryAdjacencyIndex;

    static HangingVertexBookkeeper&  getInstance();

    /**
     * @see createBoundaryVertex
     * @see createInnerVertex
     */
    tarch::la::Vector<TWO_POWER_D,int> createVertexLinkMapForNewVertex() const;
    tarch::la::Vector<TWO_POWER_D,int> createVertexLinkMapForBoundaryVertex() const;

    /**
     * Returns a real reference to an existing hanging vertex. In return, it
     * also creates a dummy hanging node if none exists yet. This is required
     * for example on level 1, where the parent hanging vertices never do exist.
     *
     * Besides the search for the right hanging node, this operation also places
     * a marker on the result that it has been used. This hanging vertex hence
     * is not garbage collected.
     *
     * @see enterCell, e.g.
     */
    tarch::la::Vector<TWO_POWER_D,int>& getAdjacencyEntriesOfVertex(
      const tarch::la::Vector<DIMENSIONS,double>&  x,
      int                                          level
    );

    bool holdsVertex(
      const tarch::la::Vector<DIMENSIONS,double>&  x,
      int                                          level
    ) const;

    bool usedVertexInThisTraversal(
      const tarch::la::Vector<DIMENSIONS,double>&  x,
      int                                          level
    ) const;

    /**
     * !!! Create a hanging node
     *
     * Forward of creating hanging nodes. If the hanging node does not
     * exist yet, the operation creates a new one with all entries set to
     * invalid. The operation always runs over all entries of the hanging
     * node. If they are unknown, the operation inherits from the coarser
     * grid. We do not only inherit in the very first iteration, but we
     * update the entries permanently. This way, we ensure that the
     * adjacency lists are updates, even if the coarser grid has not been
     * initialised completely before.
     *
     * !!! Overwrite adjacency information
     *
     * We state above that unknown entries are overwritten. That is not the
     * whole story. We also overwrite boundary data. If we have refinement
     * along the domain boundary, the creation of boundary vertices might
     * trigger new hanging nodes near the boundary. Those guys then 'inherit'
     * boundary flags from the boundary vertices though they are inside. If
     * they remain hanging, noone will ever overwrite these boundary entries
     * if the real boundary vertex's information is updated only later. So we
     * have to re-inherit these flags here as well.
     *
     * @see createHangingVertex()
     */
    tarch::la::Vector<TWO_POWER_D,int> createHangingVertex(
      const tarch::la::Vector<DIMENSIONS,double>&                  x,
      int                                                          level,
      const tarch::la::Vector<DIMENSIONS,int>&                     fineGridPositionOfVertex,
      const tarch::la::Vector<TWO_POWER_D_TIMES_TWO_POWER_D,int>&  adjacencyEntries
    );

    /**
     * If we destroy a cell, we have to invalidate all hanging vertices holding
     * that cell - otherwise those would refer to a cell not existing anymore.
     */
    void destroyCell(int cellIndex);

    template <class Vertex>
    void updateCellIndicesInMergeWithNeighbour(
      Vertex&   vertex
    );

    void beginIteration();

    /**
     * Remove all those hanging vertex entries that were not used in the previous iteration.
     */
    void endIteration();

    static bool allAdjacencyInformationIsAvailable(const tarch::la::Vector<TWO_POWER_D,int>&  arg);
    static bool allAdjacencyInformationIsAvailable(const tarch::la::Vector<THREE_POWER_D,int>&  arg);
    static bool allAdjacencyInformationIsAvailable(const tarch::la::Vector<TWO_POWER_D_TIMES_TWO_POWER_D,int>&  arg);
};


#include "multiscalelinkedcell/HangingVertexBookkeeper.cpph"


#endif
