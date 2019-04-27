// ____________________________________________________________________________
//
//  General Information:
//
//  File Name:      KShortestPathStructs.h
//  Author:         Timothy Hahn
//  Project:        KShortestPath
//
//  Description:    Declaration of structs used in the KShortestPath library
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

#ifndef _KSHORTESTPATHSTRUCTS_H_
#define _KSHORTESTPATHSTRUCTS_H_

struct kShortestPathEdges {
  size_t src_node;
  size_t dest_node;
  double edge_cost;
};

struct kShortestPathParms {
  size_t src_node;
  size_t dest_node;
  size_t k_paths;
  size_t total_nodes;
  size_t total_edges;
  kShortestPathEdges *edge_list;
};

struct kShortestPathReturn {
  size_t *pathinfo;
  double *pathcost;
  size_t *pathlen;
};

#endif
