// ____________________________________________________________________________
//
//  General Information:
//
//  File Name:      Edge.h
//  Author:         Timothy Hahn, PhD
//  Project:        raptor
//
//  Description:    The file contains the declaration of the Edge class.
//					The purpose of the Edge class is to connect the routers
//					together.
//
//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//  Revision History:
//
//  05/20/2009	v1.0	Initial Version.
//  06/02/2009	v1.02	Minor optimizations and bug fixes.
//  04/14/2019  v2.0    Reworked version based upon cmake and octave
//
// ____________________________________________________________________________

#ifndef EDGE_H
#define EDGE_H

#include <list>
#include <vector>

#include "Stats.h"

#ifdef RUN_GUI
#include "AllegroWrapper.h"
#endif

enum EdgeStatus {
	EDGE_FREE,
	EDGE_USED
};

class Edge
{
	public:
		Edge();
		Edge(int src, int dest, int spans);

		~Edge();

		inline int getSourceIndex() 
			{ return sourceIndex; };
		inline int getDestinationIndex() 
			{ return destinationIndex; };
		inline void setSourceIndex(short int s)
			{	sourceIndex = s;	};
		inline void setDestinationIndex(short int d)
			{	destinationIndex = d;	};

		inline int getNumberOfSpans() 
			{ return numberOfSpans; };
		inline void setNumberOfSpans(int s) 
			{ numberOfSpans = s; };

		inline EdgeStatus getStatus(int w)
			{ return status[w]; };
		inline int getActiveSession(int w)
			{ return activeSession[w]; };

		inline void setUsed(int session, int w)
			{ status[w] = EDGE_USED; activeSession[w] = session; };
		inline void setFree(int w)
			{ status[w] = EDGE_FREE; activeSession[w] = -1; degredation[w] = 0.0; };

		void updateUsage();
		
		void updateQMDegredation(int ci, int wavelength);
		void updateQFactorStats(int ci, int wavelength);

		std::list <void*> establishedConnections;

#ifdef RUN_GUI
		void updateGUI();
		void refreshbmps(bool useThread);
		void initializetopobmps();
#endif

		inline float getAlgorithmUsage()
			{ return algorithmUsage; };
		inline void resetAlgorithmUsage()
			{ algorithmUsage = 0.0; };
		void resetQMDegredation();

		inline float getQMDegredation()
			{ return QMDegredation; };

		inline double getPheremone()
			{ return pheremone; };
		 void resetPheremone(int ci, int spans);

		void evaporatePheremone(int ci);
		void addPheremone(int hops, int ci);

#ifdef RUN_GUI
		inline int getMaxActualUsage()
			{ return max_actual_usage; };
		inline int getUsageNums()
			{ return static_cast<int>(usageList.size()); };
		inline void addUsage(int u)
			{ usageList.push_back(u); }
		inline void setMaxUsage(int u)
			{ max_actual_usage = u; }
		void scaleEdgesTo(int spns, int px);
		void paintSpans();
		void paintUsage(int p);
		void paintUsageHistory(int x1,int y1, int x2,int y2,int t);
		void calculateAverageUsage();
		void saveData(char* file);
#endif	
		inline EdgeStats* getEdgeStats()
			{ return stats; };

		void resetEdgeStats();

		inline void insertEstablishedConnection(void* ec_void)
			{ establishedConnections.push_back(ec_void); };
		void removeEstablishedConnection(void* dcpe_void);

private:
		int sourceIndex;
		int destinationIndex;
		int numberOfSpans;

		int *activeSession;
		EdgeStatus *status;

		float algorithmUsage;
		int actualUsage;

		float QMDegredation;

		double* degredation;

#ifdef RUN_GUI
		int r1x,r1y,r2x,r2y; //coordinates of routers
		int r3x,r3y,r4x,r4y; //points to draw to for edge width (these change)
		int putX,putY;
		int thdIndx;
		float invunitX, invunitY, currX, currY, oldX, oldY, average_usage, painted_usage;//inverse unit vectors
		int painted_percent; 
		int max_actual_usage;
		vector<int> usageList;
		BITMAP* edgeBmps[14];
		BITMAP* edgeBmp;
#endif
		EdgeStats *stats;

		double pheremone;
};

#endif
