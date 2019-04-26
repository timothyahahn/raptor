// ____________________________________________________________________________
//
//  General Information:
//
//  File Name:      Router.h
//  Author:         Timothy Hahn, PhD
//  Project:        raptor
//
//  Description:    The file contains the declaration of the Router class.
//					The purpose of the Router is to simulate the optical
//					routers in the network, passing messages on their way
//					to the destination and handling them when they get to
//					the final destination.
//
//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//  Revision History:
//
//  05/20/2009	v1.0	Initial Version.
//  06/02/2009	v1.02	Minor optimizations and bug fixes.
//  04/14/2019  v2.0    Reworked version based upon cmake and octave
//
// ____________________________________________________________________________

#ifndef ROUTER_H
#define ROUTER_H

#include "Edge.h"

#ifndef NO_ALLEGRO

#include "allegro5/allegro.h"

extern int realTopoWidthPx;
extern int realTopoHeightPx;
#endif

#include <vector>

struct DP_node
{
	Edge** paths;
	bool* waveAvailability;
	unsigned int* pathLength;
	size_t* pathSpans;
	unsigned int* optimalWave;
	double* pathQuality;
	double* pathWeight;
};

class Router
{
	public:
		Router();
		~Router();

		inline void setIndex(size_t i)
			{ routerIndex = i; };
		inline size_t getIndex()
			{ return routerIndex; };

		void generateProbabilities();
		unsigned int generateDestination(double p);

#ifndef NO_ALLEGRO
		inline void incNumWorkstations()
			{ ++numWorkstations; };
		inline int getNumWorkstations()
			{ return numWorkstations; };
		void setRadius(double pct);
		inline void addTopoEdge(Edge* e)
		{	edgeList.push_back(e);	};
		void scaleEdgesTo(int spns, int px);
		inline Edge* getTopoEdgeByDest(unsigned int e)
		{	for(int x = 0; x < edgeList.size(); x++)
			{	
				if(edgeList[x]->getDestinationIndex() == e)
				{	return edgeList[x]; }	
			}
			return 0;
		};
#endif

		void addEdge(Edge* e);

		inline long long int isAdjacentTo(size_t r)
			{ return adjacencyList[r]; };

		inline Edge* getEdgeByIndex(size_t e)
			{ return edgeList[e]; };

		Edge* getEdgeByDestination(size_t r);

		void updateUsage();
		void resetUsage();

		void resetQMDegredation();
		void resetFailures();

		inline unsigned int getQualityFailures()
			{ return qualityFailures; };
		inline unsigned int getWaveFailures()
			{ return waveFailures; };

		inline void incrementQualityFailures()
			{ ++qualityFailures; };
		inline void incrementWaveFailures()
			{ ++waveFailures; };

		inline unsigned int getNumberOfEdges()
			{ return static_cast<unsigned int>(edgeList.size()); };

		void generateACOProbabilities(size_t dest);
		Edge* chooseEdge(double p);

#ifndef NO_ALLEGRO
		void refreshedgebmps(bool useThread);
		inline int getXPercent()
			{ return xpercent; };
		inline int getYPercent()
			{ return ypercent; };
		inline int getXPixel()
			{ return xPixel;	}
		inline int getYPixel()
			{ return yPixel;	}
		inline int getRadius()
		{	return radius;	};
		inline void setSelected()
		{
			if(isSelected)
				isSelected = false;
			else
				isSelected = true;
		}
		void setNumWorkstations(int w);
		void paintNumWKs();
		inline void setXPercent(int x)
		{	xpercent = x;	};
		inline void setYPercent(int y)
		{	ypercent = y;	};
		inline void setXYPixels()
		{
			xPixel = (((double)xpercent / 100.0) * realTopoWidthPx);
			yPixel = (((double)ypercent / 100.0) * realTopoHeightPx);
		}
		inline void moveXYPixels(int x, int y)
		{
			xPixel += x;
			yPixel += y;
		}
		void paintNumDests(int n);
		void paintEdgeSpans();
		void paintUsage(int p);
		void updateGUI();
		void paintProgress(int x,int y,int h, double PaintTime,double TimePerPx,int r, int g, int b);
		void saveData(char* file);

		inline void incConnAttemptsFrom()
		{	connAttemptsFromThis++;	};
		inline void setConnAttemptsFrom(unsigned int a)
		{	connAttemptsFromThis = a;	};
		inline unsigned int getConnAttemptsFrom()
		{	return connAttemptsFromThis;	};
		inline void incConnAttemptsTo()
		{	connAttemptsToThis++;	};
		inline void setConnAttemptsTo(unsigned int a)
		{	connAttemptsToThis = a;	};
		inline unsigned int getConnAttemptsTo()
		{	return connAttemptsToThis;	};
		inline void incConnSuccessesFrom()
		{	connSuccessesFromThis++;	};
		inline void setConnSuccessesFrom(unsigned int a)
		{	connSuccessesFromThis = a;	};
		inline unsigned int getConnSuccessesFrom()
		{	return connSuccessesFromThis;	};
		inline void incConnSuccessesTo()
		{	connSuccessesToThis++;	};
		inline void setConnSuccessesTo(unsigned int a)
		{	connSuccessesToThis = a;	};
		inline unsigned int getConnSuccessesTo()
		{	return connSuccessesToThis;	};
		inline void addToAvgQTo(double q)
		{	avgQTo += q;	};
		inline void addToAvgQFrom(double q)
		{	avgQFrom += q;	};

		inline const char* getName()
		{	return name;	};
		inline void setName(const char* nm)
		{	sprintf(name,nm);	};

		void selectScreen();

#endif
		DP_node* dp_node;

	private:
		size_t routerIndex;

		unsigned int qualityFailures;
		unsigned int waveFailures;

		long long int *adjacencyList;

#ifndef NO_ALLEGRO
		unsigned int connAttemptsFromThis;
		unsigned int connAttemptsToThis;
		unsigned int connSuccessesFromThis;
		unsigned int connSuccessesToThis;
		
		double avgQTo;
		double avgQFrom;

		ALLEGRO_BITMAP *routerpic;
		char name[20];
		int radius;
		int rcolor;
		int numWorkstations;
		int xpercent;
		int ypercent;
		int xPixel;
		int yPixel;
		bool isSelected;
#endif

		std::vector<Edge*> edgeList;

		double* destinationProbs;
		double* acoProbs;
};

#endif
