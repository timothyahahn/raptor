// ____________________________________________________________________________
//
//  General Information:
//
//  File Name:      QYKShortestPaths.h
//  Author:         Yan Qi
//  Project:        KShortestPath
//
//  Description:    Declaration of class(es) CQYKShortestPaths
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

#ifndef _QYKSHORTESTPATHS_H_
#define _QYKSHORTESTPATHS_H_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "QYShortestPath.h"

class CQYKShortestPaths  
{
public:
	CQYKShortestPaths(const CQYDirectedGraph& rGraph, size_t nSource, size_t nTerminal, size_t nTopk);
	virtual ~CQYKShortestPaths();

	std::vector<CQYDirectedPath*> GetTopKShortestPaths();

private: // methods

	void _Init();
	void _SearchTopKShortestPaths();

	void _DetermineCost2Target(std::vector<size_t> vertices_list, size_t deviated_node_id);
	void _RestoreEdges4CostAjustment(std::vector<size_t> vertices_list, size_t start_node_id, size_t end_node_id, bool is_deviated_node = false);
	void _UpdateWeight4CostUntilNode(size_t node_id);
	void _ReverseEdgesInGraph(CQYDirectedGraph& g);
	bool _EdgeHasBeenUsed(size_t start_node_id, size_t end_node_id);

private: // members
		
	size_t m_nTopK;
	size_t m_nSourceNodeId;
	size_t m_nTargetNodeId;
		
	const CQYDirectedGraph& m_rGraph;
	CQYDirectedGraph* m_pIntermediateGraph;
	CQYShortestPath* m_pShortestPath4IntermediateGraph;

	// variable to store the top shortest paths
	std::vector<CQYDirectedPath*> m_vTopShortestPaths;

	// a queue of candidates
	std::set<CQYDirectedPath*, CQYDirectedPath::Comparator> m_candidatePathsSet;  

	// index for node where the path derives from others
	std::map<size_t,size_t> m_pathDeviatedNodeMap;
}; 

#endif //_QYKSHORTESTPATHS_H_
