// ____________________________________________________________________________
//
//  General Information:
//
//  File Name:      ShortetsPath.cpp
//  Author:         Yan Qi
//  Project:        KShortestPath
//
//  Description:    Implementation of class(es) ShortestPath
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

#include <iostream>
#include <boost/graph/dijkstra_shortest_paths.hpp>

#include "ShortestPath.h"

/* Default Constructor
/************************************************************************/
ShortestPath::ShortestPath( const DirectedGraph& rGraph ) : m_rGraph(rGraph)
{
	m_nSourceNodeId = -1; // the source node id is 0 by default.
	_Init();
}

ShortestPath::~ShortestPath()	{}


/************************************************************************/
/* Initiate members
/************************************************************************/
void ShortestPath::_Init()
{
	size_t vertices_count = m_rGraph.GetNumberOfVertices();

	// First: edges with weights
	for (size_t i=0; i!=vertices_count; ++i)
	{
		for (size_t j=0; j!=vertices_count; ++j)
		{
			if (m_rGraph.GetWeight(i,j) != DirectedGraph::DISCONNECT)
			{
				m_vEdges.push_back(Edge_Type(i,j));
				m_vWeights.push_back(m_rGraph.GetWeight(i,j));
			}
		}
	}
}

/************************************************************************/
/* Analysis of m_vResult4Vertices and m_vResult4Distance to generate the
/* shortest path.
/************************************************************************/
DirectedPath* ShortestPath::_GetShortestPath( size_t nTargetNodeId )
{
	std::vector<size_t> vertex_list;

	// Check the input
	if (nTargetNodeId >= m_rGraph.GetNumberOfVertices() || nTargetNodeId < 0)
	{
		return new DirectedPath(-1, DirectedGraph::DISCONNECT, vertex_list);
	}

	if (m_distanceMap[nTargetNodeId] == DirectedGraph::DISCONNECT)
	{
		return new DirectedPath(-2, DirectedGraph::DISCONNECT, vertex_list);
	}

	// Determine the shortest path from the source to the terminal.
	size_t cur_vertex = nTargetNodeId;
	std::list<size_t> tmp_list;
	tmp_list.push_front(nTargetNodeId);
	do
	{
		if (m_nextNodeMap[cur_vertex] == m_nSourceNodeId)
		{
			if(cur_vertex != m_nSourceNodeId) tmp_list.push_front(m_nSourceNodeId);
			break;
		}else
		{
			cur_vertex = m_nextNodeMap[cur_vertex];
			tmp_list.push_front(cur_vertex);
		}
	} while(1);
	//
	copy(tmp_list.begin(), tmp_list.end(), back_inserter(vertex_list));

	//
	return new DirectedPath(0, m_distanceMap[nTargetNodeId], vertex_list);
}

/************************************************************************/
/* Calculate the shortest path from a source to a target.
/************************************************************************/
DirectedPath* ShortestPath::GetShortestPath( size_t nSourceNodeId, size_t nTargetNodeId )
{
	if (m_nSourceNodeId != nSourceNodeId)
	{
		m_nSourceNodeId = nSourceNodeId;
		_DijkstraShortestPathsAlg();
	}

	return _GetShortestPath(nTargetNodeId);
}

/************************************************************************/
/* Based on the input - the source of the path, create a steiner tree. (???)
/************************************************************************/
void ShortestPath::ConstructPathTree( size_t nSourceNodeId )
{
	m_nSourceNodeId = nSourceNodeId;
	_DijkstraShortestPathsAlg();
}

/************************************************************************/
/*
/************************************************************************/
void ShortestPath::_DijkstraShortestPathsAlg()
{
	size_t edges_count = m_rGraph.GetNumberOfEdges();
	size_t vertices_count = m_rGraph.GetNumberOfVertices();

	//////////////////////////////////////////////////////////////////////////
	// Initiate the boost graph
	std::vector<Vertex_Descriptor> vResult4Vertices;
	std::vector<double> vResult4Distance;
	Boost_Graph_Type g(vertices_count);
	boost::property_map<Boost_Graph_Type, boost::edge_weight_t>::type weightmap = get(boost::edge_weight, g);
	//
	for (std::size_t j = 0; j < static_cast<size_t>(edges_count); ++j)
	{
		Edge_Descriptor e; bool inserted;
		tie(e, inserted) = add_edge(m_vEdges[j].first, m_vEdges[j].second, g);
		weightmap[e] = m_vWeights[j];
	}

	// about the vertices in the boost graph
	vResult4Vertices.resize(num_vertices(g));
	vResult4Distance.resize(num_vertices(g));
	Vertex_Descriptor s = vertex(m_nSourceNodeId, g);

	// run the algorithm
	// VC++ has trouble with the named parameters mechanism
	boost::property_map<Boost_Graph_Type, boost::vertex_index_t>::type indexmap = get(boost::vertex_index, g);
	dijkstra_shortest_paths(g, s, &vResult4Vertices[0], &vResult4Distance[0], weightmap, indexmap,
		std::less<double>(), boost::closed_plus<double>(),
		DirectedGraph::DISCONNECT, 0, boost::default_dijkstra_visitor());

	//////////////////////////////////////////////////////////////////////////
	// Set the results
	for (size_t i = 0; i < vResult4Vertices.size(); ++i)
	{
		m_distanceMap[i] = vResult4Distance[i];
		m_nextNodeMap[i] = vResult4Vertices[i];
	}
}
