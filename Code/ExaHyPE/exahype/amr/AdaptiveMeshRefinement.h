/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released under the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 *
 *
 * TODO(Dominic): We currently just collect methods here.
 * These methods are not or only partially used at the moment.
 **/

#ifndef ADAPTIVEMESHREFINEMENT_H_
#define ADAPTIVEMESHREFINEMENT_H_

#include "peano/utils/Globals.h"

#include "exahype/Cell.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

namespace exahype {

namespace amr {
  /**
   * Returns  AugmentationControl::::NextToCell if the cell has a neighbour of
   * type exahype::records::ADERDGCellDescription:Cell,
   * exahype::solvers::Solver::NextToAncestor if the cell has a neighbour of
   * type exahype::records::ADERDGCellDescription:Ancestor or
   * exahype::records::ADERDGCellDescription::EmptyAncestor.
   * Is both the case, this function returns
   * AugmentationControl::::NextToCellOrAncestor.
   * If none of the previous is the case, this function returns
   * AugmentationControl::Default.
   */
  template<class CellDescription,class CellDescriptionsHeap>
  exahype::solvers::Solver::AugmentationControl
  augmentationCriterion(
      const int solverNumber,
      const typename CellDescription::Type type,
      const int level,
      const tarch::la::Vector<THREE_POWER_D, int>&
      neighbourCellDescriptionsIndices) {
    // left,right,front,back,(bottom,top)
#if DIMENSIONS == 2
    constexpr int neighbourPositions[4] = {3, 5, 1, 7};
#else
    constexpr int neighbourPositions[6] = {12, 14, 10, 16, 4, 22};
#endif
    bool nextToAncestor = false;
    bool nextToCell = false;

    for (int i = 0; i < DIMENSIONS_TIMES_TWO; i++) {
      const int neighbourCellDescriptionIndex =
          neighbourCellDescriptionsIndices[neighbourPositions[i]];
      if (CellDescriptionsHeap::getInstance().isValidIndex(neighbourCellDescriptionIndex)) {
        for (CellDescription& pNeighbour : CellDescriptionsHeap::getInstance().getData(
            neighbourCellDescriptionIndex)) {
          if (pNeighbour.getSolverNumber() == solverNumber &&
              pNeighbour.getLevel() == level) {
            switch (pNeighbour.getType()) {
              case CellDescription::Ancestor:
              case CellDescription::EmptyAncestor:
                nextToAncestor = true;
                break;
              case CellDescription::Cell:
                nextToCell = true;
                break;
              default:
                break;
            }
          }
        }
      }
    }

    // NOTE: The order below is important.
    if (nextToCell && nextToAncestor) {
      return exahype::solvers::Solver::AugmentationControl::NextToCellAndAncestor;
    }
    if (nextToAncestor) {
      return exahype::solvers::Solver::AugmentationControl::NextToAncestor;
    }
    if (nextToCell) {
      return exahype::solvers::Solver::AugmentationControl::NextToCell;
    }
    // Erase otherwise.
    return exahype::solvers::Solver::AugmentationControl::Default;
  }

  /*
   * Change the erasing request to a change to descendant request if the coarse grid Cell
   * has children (of type Descendant).
   * Rationale: We cannot directly erase a Cell that has children (of type Descendant).
   *
   * Further, reset the deaugmenting request if a coarse grid Descendant has children
   * (of type Descendant). Rationale:
   * We cannot erase a coarse grid cell that has children (of type Descendant)
   * before erasing the children.
   *
   * \note This method should be called when we enter a fine grid cell,
   * coarseGridCellDescription is a cell description associated with
   * the coarse grid cell.
   *
   * \note A more sophisticated procedure has to performed for the refinement event
   * AugmentationRequested. We need to use the taversal's descend event to handle
   * this event. We thus do not rely on fineGridCell.isRefined() in the previous enterCell event
   * to check if we need to reset the deaugmenting request.
   */
  template<class CellDescription>
  void resetErasingOrDeaugmentingRequestIfParent(
      CellDescription& coarseGridCellDescription) {
    switch (coarseGridCellDescription.getRefinementEvent()) {
      case CellDescription::DeaugmentingRequested:
        assertion1(coarseGridCellDescription.getType()==CellDescription::EmptyDescendant ||
                   coarseGridCellDescription.getType()==CellDescription::Descendant,coarseGridCellDescription.toString());
        coarseGridCellDescription.setRefinementEvent(CellDescription::None);
        break;
      case CellDescription::ErasingRequested:
        assertion1(coarseGridCellDescription.getType()==CellDescription::Cell,
                   coarseGridCellDescription.toString());
        coarseGridCellDescription.setRefinementEvent(CellDescription::ChangeToDescendantRequested);
        break;
      default:
        break;
    }
  }

  /**
   * Determine the position of a Cell or Ancestor with respect
   * to a parent of type Ancestor.
   *
   * This method is required for the face data restriction, the
   * volume data restriction, and the FV volume data restriction.
   *
   * @todo:16/04/09:Dominic Etienne Charrier: I am not sure if this the
   * right file/class to hold this functionality.
   */
  template <class CellDescription,class CellDescriptionHeap>
  exahype::solvers::Solver::SubcellPosition
  computeSubcellPositionOfCellOrAncestor(
      const CellDescription& pChild) {
    assertion1(pChild.getType()==CellDescription::Cell ||
               pChild.getType()==CellDescription::Ancestor,pChild.toString());

    exahype::Cell::SubcellPosition subcellPosition;
    // Initialisation.
    subcellPosition.parentIndex = pChild.getParentIndex();
    CellDescription* pParent = 0;
    for (auto& p : CellDescriptionHeap::getInstance().getData(
        subcellPosition.parentIndex)) {  // Loop over cell descriptions
      if (p.getSolverNumber() == pChild.getSolverNumber()) {
        pParent = &p;
        assertion(CellDescriptionHeap::getInstance().isValidIndex(
            subcellPosition.parentIndex));
      }
    }

    if (pParent != 0) {
      // Iterative determining of the top most parent that might hold data.
      while (pParent->getType() == CellDescription::EmptyAncestor &&
          CellDescriptionHeap::getInstance().isValidIndex(pParent->getParentIndex())) {
        const int currentParentIndex =
            pParent->getParentIndex();
        // Value must be fixed. We update pParent within the loop.
        for (auto& p : CellDescriptionHeap::getInstance().getData(
            currentParentIndex)) {  // Loop over cell descriptions
          if (p.getSolverNumber() == pChild.getSolverNumber()) {
            subcellPosition.parentIndex = pParent->getParentIndex();
            pParent = &p;
          }
        }
      }
      assertion(pParent->getType()==CellDescription::Ancestor ||
                pParent->getType()==CellDescription::EmptyAncestor);

      // Compute subcell index.
      for (int xi = 0; xi < DIMENSIONS; ++xi) {
        assertion((pChild.getOffset(xi) - pParent->getOffset(xi)) >= 0);
        subcellPosition.subcellIndex[xi] = tarch::la::round(
            (pChild.getOffset(xi) - pParent->getOffset(xi))/pChild.getSize(xi));
      }
    }

    return subcellPosition;
  }

  /**
   * Determine the position of a Descendant with respect
   * to a  Cell or Descendant that contains data, i.e.,
   * has at least one neighbour that is a real cell.
   *
   * This method is required for the face data prolongation, the
   * volume data prolongation, and the FV volume data prolongation.
   *
   * @todo:16/04/09:Dominic Etienne Charrier: I am not sure if this the
   * right file/class to hold this functionality.
   */
  template <class CellDescription,class CellDescriptionHeap>
  exahype::solvers::Solver::SubcellPosition
  computeSubcellPositionOfDescendant(
      const CellDescription& pChild) {
    assertion1(pChild.getType()==CellDescription::Descendant,pChild.toString());

    exahype::Cell::SubcellPosition subcellPosition;

    // Initialisation.
    assertion1(CellDescriptionHeap::getInstance().isValidIndex(
        pChild.getParentIndex()),
        pChild.getParentIndex());
    subcellPosition.parentIndex = pChild.getParentIndex();
    CellDescription* pParent = 0;

    for (auto& p : CellDescriptionHeap::getInstance().getData(
        pChild.getParentIndex())) {
      if (p.getSolverNumber() == pChild.getSolverNumber()) {
        pParent = &p;
      }
    }

    if (pParent != 0) {
      // recursion
      while (pParent->getType() == CellDescription::EmptyDescendant) {
        const int currentParentIndex =
            pParent->getParentIndex();
        // Value must be fixed. We update pParent within the loop.
        assertion1(CellDescriptionHeap::getInstance().isValidIndex(
            currentParentIndex),
            currentParentIndex);
        for (auto& p : CellDescriptionHeap::getInstance().getData(currentParentIndex)) {
          if (p.getSolverNumber() == pChild.getSolverNumber()) {
            subcellPosition.parentIndex = pParent->getParentIndex();
            pParent = &p;
          }
        }
      }

      assertion(pParent->getType() == CellDescription::Descendant ||
                pParent->getType() == CellDescription::Cell);

      // Ccompute subcell index.
      for (int xi = 0; xi < DIMENSIONS; ++xi) {
        assertion((pChild.getOffset(xi) - pParent->getOffset(xi)) >= 0);
        subcellPosition.subcellIndex[xi] = tarch::la::round(
            (pChild.getOffset(xi) - pParent->getOffset(xi))/pChild.getSize(0));
      }
    } else {
      std::cerr << "computeSubcellPositionOfDescendant: parent of "
          "descendant could not be found!"
          << std::endl;
      exit(EXIT_FAILURE);
    }

    return subcellPosition;
  }

}  // namespace amr
}  // namespace exahype



#endif /* ADAPTIVEMESHREFINEMENT_H_ */
