// ____________________________________________________________________________
//
//  General Information:
//
//  File Name:      ShortestPath.h
//  Author:         Yan Qi
//  Project:        KShortestPath
//
//  Description:    Declaration of class ShortestPath, which implements
//  Dijkstra algorithm for the shortest path in the directed graph.
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

#ifndef _QYSHORTETSPATH_H_
#define _QYSHORTETSPATH_H_

#include <map>
#include <vector>

#include "ConfigCenter.h"
#include "DirectedGraph.h"
#include "DirectedPath.h"

class ShortestPath {
 public:
  ShortestPath(const DirectedGraph& rGraph);
  virtual ~ShortestPath();

  DirectedPath* GetShortestPath(size_t nSourceNodeId, size_t nTargetNodeId);
  void ConstructPathTree(size_t nSourceNodeId);

  double GetDistance(size_t i) {
    return m_distanceMap.count(i) ? m_distanceMap[i]
                                  : DirectedGraph::DISCONNECT;
  }
  void SetDistance(size_t i, double new_value) { m_distanceMap[i] = new_value; }

  size_t GetNextNodeId(size_t i) {
    return m_distanceMap.count(i) ? m_nextNodeMap[i] : DirectedGraph::DEADEND;
  }
  void SetNextNodeId(size_t i, size_t val) { m_nextNodeMap[i] = val; }

 private:  // methods
  void _Init();
  void _DijkstraShortestPathsAlg();
  DirectedPath* _GetShortestPath(size_t nTargetNodeId);

 private:  // members
  std::vector<ConfigCenter::SizeT_Pair> m_vEdges;
  std::vector<double> m_vWeights;

  std::map<size_t, double> m_distanceMap;
  std::map<size_t, size_t> m_nextNodeMap;

  size_t m_nSourceNodeId;
  const DirectedGraph& m_rGraph;
};

#endif  //_QYSHORTETSPATH_H_
