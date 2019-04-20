// ____________________________________________________________________________
//
//  General Information:
//
//  File Name:      CQYShortestPath.h
//  Author:         Timothy Hahn
//  Project:        KShortestPath
//
//  Description:    Declaration of class CQYShortestPath, which implements 
//  Dijkstra algorithm for the shortest path in the directed graph.
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

#ifndef _QYINCLUDE_H_
#define _QYINCLUDE_H_

struct kShortestPathEdges {
	unsigned short int src_node;
	unsigned short int dest_node;
	float edge_cost;
};

struct kShortestPathParms {
	unsigned short int src_node;
	unsigned short int dest_node;
	unsigned short int k_paths;
	unsigned short int total_nodes;
	unsigned short int total_edges;
	kShortestPathEdges *edge_list;	
};

struct kShortestPathReturn {
	unsigned short int *pathinfo;
	float *pathcost;
	unsigned short int *pathlen;
};

#endif
