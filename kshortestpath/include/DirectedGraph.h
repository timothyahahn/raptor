// ____________________________________________________________________________
//
//  General Information:
//
//  File Name:      DirectedGraph.h
//  Author:         Yan Qi
//  Project:        KShortestPath
//
//  Description:    Declaration of class(es) DirectedGraph
//
//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//  Revision History:
//
//  11/21/2006   Yan   Initial Version
//  01/11/2007   Yan   Modified Version: correct the way to calculate the number
//  of edges in the graph. 04/19/2019   Hahn  Modern OS/Compiler Changes
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

#ifndef _DIRECTEDGRAPH_H_
#define _DIRECTEDGRAPH_H_

#include <string>

#include "ConfigCenter.h"
#include "KShortestPathStructs.h"

// ____________________________________________________________________________
//
// Class:       DirectedGraph
//
// Purpose:     DirectedGraph defines the directed graph with a list of
//              directed edges, associated with its weight.
//
// Notes:		Two ways to construct a graph:
//				1. Assign the path of the file recording the
//graph
//				2. Transfer an existing object.
//
// See Also:
//
// ____________________________________________________________________________
class DirectedGraph {
 public:
  const static double DISCONNECT;
  const static size_t DEADEND;

  DirectedGraph();
  DirectedGraph(kShortestPathParms params);
  virtual ~DirectedGraph();

  DirectedGraph(const DirectedGraph& rGraph);
  DirectedGraph& operator=(const DirectedGraph& rGraph);

  void RemoveEdge(size_t i, size_t j);
  void AddEdge(size_t i, size_t j, double weight);

  // Getter and setter
  size_t GetNumberOfVertices() const { return m_nNumberOfVertices; }
  void SetNumberOfVertices(size_t val) { m_nNumberOfVertices = val; }

  size_t GetNumberOfEdges() const { return m_nNumberOfEdges; }
  void SetNumberOfEdges(size_t val) { m_nNumberOfEdges = val; }

  double GetMaxWeight() const { return m_dMaxWeight; }
  void SetMaxWeight(double val) { m_dMaxWeight = val; }

  double GetMinWeight() const { return m_dMinWeight; }
  void SetMinWeight(double val) { m_dMinWeight = val; }

  double GetWeight(size_t i, size_t j) const {
    return m_pDirectedEdges->count(ConfigCenter::SizeT_Pair(i, j))
               ? (*m_pDirectedEdges)[ConfigCenter::SizeT_Pair(i, j)]
               : DISCONNECT;
  }
  void SetWeight(size_t i, size_t j, double val) {
    (*m_pDirectedEdges)[ConfigCenter::SizeT_Pair(i, j)] = val;
  }

 private:
  ConfigCenter::SizeT_Pair_Double_Map* m_pDirectedEdges;

  size_t m_nNumberOfVertices;
  size_t m_nNumberOfEdges;

  double m_dMaxWeight;
  double m_dMinWeight;

  void _Init();
};

#endif  //_DIRECTEDGRAPH_H_
