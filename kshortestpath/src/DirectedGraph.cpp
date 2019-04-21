// ____________________________________________________________________________
//
//  General Information:
//
//  File Name:      DirectedGraph.cpp
//  Author:         Yan Qi
//  Project:        KShortestPath
//
//  Description:    Implementation of class(es) DirectedGraph
//
//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//  Revision History:
//
//  11/21/2006   Yan   Initial Version
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

#include <limits>
#include <fstream>
#include <iostream>

#include "DirectedGraph.h"

const size_t DirectedGraph::DEADEND = -1;
const double DirectedGraph::DISCONNECT = (std::numeric_limits<double>::max)();

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
DirectedGraph::DirectedGraph():
	m_nNumberOfVertices(0), m_nNumberOfEdges(0), m_dMaxWeight(0.0), m_dMinWeight(DISCONNECT),
	m_pDirectedEdges(new ConfigCenter::SizeT_Pair_Double_Map())
{
}

DirectedGraph::DirectedGraph( kShortestPathParms params ) :
	m_nNumberOfVertices(0), m_nNumberOfEdges(0), m_dMaxWeight(0.0), m_dMinWeight(DISCONNECT),
	m_pDirectedEdges(new ConfigCenter::SizeT_Pair_Double_Map())
{
	m_nNumberOfVertices = params.total_nodes;
		
	for(size_t a = 0; a < params.total_edges; ++a)
	{
		size_t i = params.edge_list[a].src_node;
		size_t j = params.edge_list[a].dest_node;
		double w = params.edge_list[a].edge_cost;

		m_pDirectedEdges->insert(std::pair<std::pair<size_t, size_t>, double>(std::pair<size_t, size_t>(i,j), w));
			
		++m_nNumberOfEdges;
			
		if (w > m_dMaxWeight)
		{
			m_dMaxWeight = w;
		}
			
		if (w < m_dMinWeight)
		{
			m_dMinWeight = w;
		}
	}	
		
	m_nNumberOfEdges = m_pDirectedEdges->size();
}

DirectedGraph::DirectedGraph( const DirectedGraph& rGraph ) :
	m_nNumberOfVertices(0), m_nNumberOfEdges(0), m_dMaxWeight(0.0), m_dMinWeight(DISCONNECT),
	m_pDirectedEdges(new ConfigCenter::SizeT_Pair_Double_Map())
{
	*this = rGraph;	
}

DirectedGraph& DirectedGraph::operator=( const DirectedGraph& rGraph )
{
	m_nNumberOfVertices = rGraph.m_nNumberOfVertices;
	m_nNumberOfEdges = rGraph.m_nNumberOfEdges;

	m_pDirectedEdges = new ConfigCenter::SizeT_Pair_Double_Map(*(rGraph.m_pDirectedEdges));
		
	return *this;
}

DirectedGraph::~DirectedGraph()
{
	if (m_pDirectedEdges != nullptr)
	{
		delete m_pDirectedEdges;
	}
}

void DirectedGraph::_Init()
{
	m_nNumberOfEdges = 0;
	m_dMaxWeight = 0;
	m_dMinWeight = DISCONNECT;
	m_pDirectedEdges = new ConfigCenter::SizeT_Pair_Double_Map();
}

void DirectedGraph::RemoveEdge( size_t i, size_t j )
{
	ConfigCenter::SizeT_Pair_Double_Map_Iterator pos = m_pDirectedEdges->find(ConfigCenter::SizeT_Pair(i,j));
	if (pos != m_pDirectedEdges->end())
	{
		m_pDirectedEdges->erase(pos);
	}
}

void DirectedGraph::AddEdge( size_t i, size_t j, double weight )
{
	(*m_pDirectedEdges)[ConfigCenter::SizeT_Pair(i,j)] = weight;
}
