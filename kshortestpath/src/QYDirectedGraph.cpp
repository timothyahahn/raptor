// ____________________________________________________________________________
//
//  General Information:
//
//  File Name:      QYDirectedGraph.cpp
//  Author:         Yan Qi
//  Project:        KShortestPath
//
//  Description:    Implementation of class(es) CQYDirectedGraph
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

#ifdef _MSC_VER
#pragma warning(disable: 4786)
#endif // _MSC_VER 

#include <limits>
#include <fstream>
#include <iostream>
#include "QYDirectedGraph.h"

const int CQYDirectedGraph::DEADEND = -1;
const double CQYDirectedGraph::DISCONNECT = (std::numeric_limits<double>::max)();

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
CQYDirectedGraph::CQYDirectedGraph():
	m_nNumberOfVertices(0), m_nNumberOfEdges(0), m_dMaxWeight(0.0), m_dMinWeight(DISCONNECT),
	m_pDirectedEdges(new CQYConfigCenter::IntPair_Double_Map())
{
}

CQYDirectedGraph::CQYDirectedGraph( kShortestPathParms params ) :
	m_nNumberOfVertices(0), m_nNumberOfEdges(0), m_dMaxWeight(0.0), m_dMinWeight(DISCONNECT),
	m_pDirectedEdges(new CQYConfigCenter::IntPair_Double_Map())
{
	m_nNumberOfVertices = params.total_nodes;
		
	for(unsigned int a = 0; a < params.total_edges; ++a)
	{
		int i = params.edge_list[a].src_node;
		int j = params.edge_list[a].dest_node;
		double w = params.edge_list[a].edge_cost;

		m_pDirectedEdges->insert(std::pair<std::pair<int, int>, double>(std::pair<int, int>(i,j), w));
			
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

CQYDirectedGraph::CQYDirectedGraph( const CQYDirectedGraph& rGraph )
{
	*this = rGraph;	
}

CQYDirectedGraph& CQYDirectedGraph::operator=( const CQYDirectedGraph& rGraph )
{
	m_nNumberOfVertices = rGraph.m_nNumberOfVertices;
	m_nNumberOfEdges = rGraph.m_nNumberOfEdges;

	m_pDirectedEdges = new CQYConfigCenter::IntPair_Double_Map(*(rGraph.m_pDirectedEdges));
		
	return *this;
}

CQYDirectedGraph::~CQYDirectedGraph()
{
	if (m_pDirectedEdges != nullptr)
	{
		delete m_pDirectedEdges;
	}
}

void CQYDirectedGraph::_Init()
{
	m_nNumberOfEdges = 0;
	m_dMaxWeight = 0;
	m_dMinWeight = DISCONNECT;
	m_pDirectedEdges = new CQYConfigCenter::IntPair_Double_Map();
}

void CQYDirectedGraph::RemoveEdge( int i, int j )
{
	CQYConfigCenter::IntPair_Double_Map_Iterator pos = m_pDirectedEdges->find(std::pair<int, int>(i,j));
	if (pos != m_pDirectedEdges->end())
	{
		m_pDirectedEdges->erase(pos);
	}
}

void CQYDirectedGraph::AddEdge( int i, int j, double weight )
{
	(*m_pDirectedEdges)[std::pair<int, int>(i,j)] = weight;	
}
