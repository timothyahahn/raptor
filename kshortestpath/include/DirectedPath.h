// ____________________________________________________________________________
//
//  General Information:
//
//  File Name:      QYDirectedPath.h
//  Author:         Yan Qi
//  Project:        KShortestPath
//
//  Description:    Declaration of class(es) DirectedPath
//
//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//  Revision History:
//
//  11/21/2006   Yan   Initial Version
//  01/11/2008   Yan   Modified Version
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

#ifndef _QYDIRECTEDPATH_H_
#define _QYDIRECTEDPATH_H_

#include <vector>
#include <fstream>
#include <algorithm>
#include <iterator>

class DirectedPath
{
public:
	DirectedPath()
	: m_nId(0), m_nLength(0), m_dCost(0), m_nSourceNodeId(0), m_nTerminalNodeId(0)
	{
	}

	DirectedPath(size_t pId, double pCost, const std::vector<size_t>& pVertexList)
		:m_nId(pId), m_nLength(0), m_dCost(pCost), m_nSourceNodeId(0), m_nTerminalNodeId(0)
	{
		m_vVertexList.assign(pVertexList.begin(), pVertexList.end());
	}
		
	virtual ~DirectedPath(){};
		
	// Getter and Setter
	size_t GetId() const { return m_nId; }
	void SetId(size_t val) { m_nId = val; }
		
	double GetCost() const { return m_dCost; }
	void SetCost(double val) { m_dCost = val; }
		
	size_t GetLength() const { return m_vVertexList.size(); }
		
	std::vector<size_t> GetVertexList() const { return m_vVertexList; }
	void SetVertexList(std::vector<size_t> val) { m_vVertexList = val; }
		
	size_t GetSourceNodeId() const { return m_nSourceNodeId; }
	void SetSourceNodeId(size_t val) { m_nSourceNodeId = val; }
		
	size_t GetTerminalNodeId() const { return m_nTerminalNodeId; }
	void SetTerminalNodeId(size_t val) { m_nTerminalNodeId = val; }
		
	// display the content
	void PrintOut(std::ostream& out_stream) const
	{
		out_stream << "Cost: " << m_dCost << " Length: " << m_vVertexList.size() << std::endl;
		std::copy(m_vVertexList.rbegin(), m_vVertexList.rend(), std::ostream_iterator<int>(out_stream, " "));
		out_stream << std::endl <<  "*********************************************" << std::endl;	
	}


private: // members
	size_t m_nId;
	size_t m_nLength;
	double m_dCost;  
	std::vector<size_t> m_vVertexList;
		
	// intermediate variables
	size_t m_nSourceNodeId;
	size_t m_nTerminalNodeId;

public:
	//// Comparator for paths: the smaller path has less cost.
	class Comparator 
	{
	public:
		// Lesson: the condition must be checked completely!!!
		bool operator() (const DirectedPath& s1, const DirectedPath& s2) const 
		{
			if (s1.GetCost() == s2.GetCost())
			{
				if (s1.GetLength() == s2.GetLength())
				{
					return s1.GetId() < s2.GetId();
				}
				else
				{
					return s1.GetLength() < s2.GetLength();
				}
			}
			else
			{
				return s1.GetCost() < s2.GetCost();
			}
		}

		//
		bool operator() (const DirectedPath* s1, const DirectedPath* s2) const 
		{
			return operator()(*s1, *s2);
		}
	}; 
};

#endif //_QYDIRECTEDPATH_H_
