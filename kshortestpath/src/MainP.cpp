// ____________________________________________________________________________
//
//  General Information:
//
//  File Name:      DirectedGraph.cpp
//  Author:         Timothy Hahn
//  Project:        KShortestPath
//
//  Description:    Implementation of class(es) DirectedGraph
//
//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//  Revision History:
//
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

#include "ShortestPath.h"
#include "KShortestPaths.h"

#ifdef __GNUC__ 
	extern "C" void calc_k_shortest_paths(const kShortestPathParms &params, kShortestPathReturn* retVal);
#else
	extern "C" __declspec(dllexport) void calc_k_shortest_paths(const kShortestPathParms &params, kShortestPathReturn* retVal);
#endif

void copyResults(std::vector<DirectedPath*> &topK_shortest_paths, const kShortestPathParms &params, kShortestPathReturn* retVal, bool dijkstra);

void calc_k_shortest_paths(const kShortestPathParms &params, kShortestPathReturn* retVal)
{
	DirectedGraph dg(params);

	if(params.k_paths == 1)
	{
		std::vector<DirectedPath*> topK_shortest_paths;

		ShortestPath sp(dg);
		
		topK_shortest_paths.push_back(sp.GetShortestPath(params.src_node, params.dest_node));

		copyResults(topK_shortest_paths, params, retVal, true);

		if(topK_shortest_paths.size() > 0)
			delete topK_shortest_paths[0];
	}
	else
	{
		KShortestPaths ksp(dg, params.src_node, params.dest_node, params.k_paths);

		std::vector<DirectedPath*> topK_shortest_paths = ksp.GetTopKShortestPaths();

		copyResults(topK_shortest_paths, params, retVal, false);
	}

	return;
}

void copyResults(std::vector<DirectedPath*> &topK_shortest_paths, const kShortestPathParms &params, kShortestPathReturn* retVal, bool dijkstra)
{
	for(std::vector<DirectedPath*>::iterator iter = topK_shortest_paths.begin(); iter != topK_shortest_paths.end(); ++iter)
	{
		if((*iter)->GetLength() <= 0 || (*iter)->GetLength() >= params.total_nodes)
		{
			if(dijkstra == true)
				delete *iter;

			topK_shortest_paths.erase(iter);
			
			if(topK_shortest_paths.size() == 0)
				break;
			else
				iter = topK_shortest_paths.begin();
		}
	}

	for(size_t a = 0; a < topK_shortest_paths.size(); ++a)
	{
		retVal->pathcost[a] = double(topK_shortest_paths[a]->GetCost());
		retVal->pathlen[a] = topK_shortest_paths[a]->GetLength();

		for(size_t b = 0; b < retVal->pathlen[a]; ++b)
		{
			retVal->pathinfo[a * (params.total_nodes - 1) + b] = topK_shortest_paths[a]->GetVertexList()[b];
		}
	}

	for(size_t d = topK_shortest_paths.size(); d < params.k_paths; ++d)
	{
		retVal->pathcost[d] = std::numeric_limits<double>::infinity();
		retVal->pathlen[d] = std::numeric_limits<size_t>::infinity();
	}
}
