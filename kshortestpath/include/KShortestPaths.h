// ____________________________________________________________________________
//
//  General Information:
//
//  File Name:      KShortestPaths.h
//  Author:         Yan Qi
//  Project:        KShortestPath
//
//  Description:    Declaration of class(es) KShortestPaths
//
//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//  Revision History:
//
//  11/23/2006   Yan   Initial Version
//  04/19/2019   Hahn  Modern OS/Compiler Changes
//
//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//  Copyright Notice:
//
//  Copyright (c) 2006 Your Company Inc.
//
//  Warning: This computer program is protected by copyright law and
//  international treaties.  Unauthorized reproduction or distribution
//  of this program, or any portion of it, may result in severe civil and
//  criminal penalties, and will be prosecuted to the maximum extent
//  possible under the law.
//
// ____________________________________________________________________________

#ifndef _KShortestPaths_H_
#define _KShortestPaths_H_

#include <set>

#include "DirectedGraph.h"
#include "DirectedPath.h"
#include "ShortestPath.h"

class KShortestPaths {
 public:
  KShortestPaths(const DirectedGraph& rGraph, size_t nSource, size_t nTerminal,
                 size_t nTopk);
  virtual ~KShortestPaths();

  std::vector<DirectedPath*> GetTopKShortestPaths();

 private:  // methods
  void _SearchTopKShortestPaths();

  void _DetermineCost2Target(std::vector<size_t> vertices_list,
                             size_t deviated_node_id);
  void _RestoreEdges4CostAjustment(std::vector<size_t> vertices_list,
                                   size_t start_node_id, size_t end_node_id,
                                   bool is_deviated_node = false);
  void _UpdateWeight4CostUntilNode(size_t node_id);
  void _ReverseEdgesInGraph(DirectedGraph& g);
  bool _EdgeHasBeenUsed(size_t start_node_id, size_t end_node_id);

 private:  // members
  size_t m_nTopK;
  size_t m_nSourceNodeId;
  size_t m_nTargetNodeId;

  const DirectedGraph& m_rGraph;
  DirectedGraph* m_pIntermediateGraph;
  ShortestPath* m_pShortestPath4IntermediateGraph;

  // variable to store the top shortest paths
  std::vector<DirectedPath*> m_vTopKShortestPaths;

  // a queue of candidates
  std::set<DirectedPath*, DirectedPath::Comparator> m_candidatePathsSet;

  // index for node where the path derives from others
  std::map<size_t, size_t> m_pathDeviatedNodeMap;
};

#endif  //_KShortestPaths_H_
