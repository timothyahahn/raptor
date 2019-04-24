// ____________________________________________________________________________
//
//  General Information:
//
//  File Name:      Thread.h
//  Author:         Timothy Hahn, PhD
//  Project:        raptor
//
//  Description:    The file contains the declaration of the thread class.
//					The purpose of the controller is to manage the simulation and
//					its associated events. This includes controlling the simulation
//					time, managing events, controlling the active/inactive work-
//					stations, and collecting/saving/printing the results of the
//					simulation.
//
//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//  Revision History:
//
//  05/20/2009	v1.0	Initial Version.
//  06/02/2009	v1.02	Minor optimizations and bug fixes.
//  04/14/2019  v2.0    Reworked version based upon cmake and octave
//
// ____________________________________________________________________________

#ifndef CONTROLLER_H
#define CONTROLLER_H

#include <cmath>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <limits>
#include <queue>
#include <random>
#include <string>
#include <vector>

#include "AlgorithmParameters.h"
#include "EstablishedConnections.h"
#include "ErrorCodes.h"
#include "EventQueue.h"
#include "MessageLogger.h"
#include "QualityParameters.h"
#include "ResourceManager.h"
#include "Router.h"
#include "Stats.h"
#include "Workstation.h"

class Thread
{
	public:
		Thread();
		Thread(int ci, int argc, const char* argv[], bool isLPS);
		~Thread();

		inline Router* getRouterAt(size_t i)
			{ return routers[i]; };
		inline Workstation* getWorkstationAt(size_t i)
			{ return workstations[i]; };

		inline void addRouter(Router* r)
			{ routers.push_back(r); };
		inline void addWorkstation(Workstation* w)
			{ workstations.push_back(w); };

		inline double getGlobalTime()
			{ return globalTime; };
		inline void setGlobalTime(double t)
			{ globalTime = t; };
		
		inline GlobalStats& getGlobalStats()
			{ return stats; };

		inline void setControllerIndex(unsigned int ci)
			{ controllerIndex = ci; };

		void initThread(AlgorithmToRun* alg);

		int runThread(AlgorithmToRun* alg);

		void initPriorityQueue(unsigned int w);
		void initResourceManager();

		inline QualityParameters getQualityParams()
			{ return qualityParams; };

		inline void recordEvent(const std::string& s, bool print, unsigned int ci)
		{
			if(isLoadPrevious == true)
				return;
			else if(controllerIndex == 0)
				logger->recordEvent(s,print,ci);
			else
				exit(ERROR_RECORD_EVENT);
		};

		inline  void flushLog(bool print)
		{
			if (isLoadPrevious == true)
				return;
			else if (controllerIndex == 0)
				logger->flushLog(print);
			else
				exit(ERROR_NO_FLUSH);
		}

		inline size_t getNumberOfRouters()
			{ return numberOfRouters; };
		inline size_t getNumberOfWorkstations()
			{ return numberOfWorkstations; };
		inline size_t getNumberOfEdges()
			{ return numberOfEdges; };
		inline size_t getNumberOfWavelengths()
			{ return numOfWavelengths; };
		inline void setNumberOfWavelengths(size_t n)
			{	numOfWavelengths = n;	};
		inline size_t getRandomSeed()
			{ return randomSeed; };

		inline ResourceManager* getResourceManager()
			{ return rm; };

		inline RoutingAlgorithm getCurrentRoutingAlgorithm()
			{ return CurrentRoutingAlgorithm; };
		inline WavelengthAlgorithm getCurrentWavelengthAlgorithm()
			{ return CurrentWavelengthAlgorithm; };
		inline ProbeStyle getCurrentProbeStyle()
			{ return CurrentProbeStyle; };
		inline bool getCurrentQualityAware()
			{ return CurrentQualityAware; };
		inline unsigned int getCurrentActiveWorkstations()
			{ return CurrentActiveWorkstations; };
		inline void setCurrentActiveWorkstations(unsigned int w)
			{ CurrentActiveWorkstations = w;	}

		double getBeta();

		inline double getMinDuration()
			{ return minDuration; };
		inline size_t getMaxSpans()
			{ return maxSpans; };

		void setMinDuration(size_t spans);
		void setQFactorMin(size_t spans);

		inline const std::string getRoutingAlgorithmName(unsigned int a)
			{ return RoutingAlgorithmNames[a]; };
		inline const std::string getWavelengthAlgorithmName(unsigned int w)
			{ return WavelengthAlgorithmNames[w]; };
		inline const std::string getProbeStyleName(unsigned int p)
			{ return ProbeStyleNames[p]; };

#ifdef RUN_GUI
		void setTerminate()
			{ terminateProgram = true; };

		void saveThread(char* dir);

		inline void setMaxMaxUsage(double m)
			{ maxMaxUsage = m; };
		inline double getMaxMaxUsage()
			{ return maxMaxUsage; };

		inline void setPaintDir()
		{	if( paintDir )
			{	paintDir = false;	}
			else
			{	paintDir = true;	}
		};

		inline bool getPaintDir()
			{	return paintDir; };

		inline void setPaintSpans()//I don't think this should be here. Do something similar
		{	if( paintSpans )		//in GUI.cpp, instead
			{	paintSpans = false;	}
			else
			{	paintSpans = true;	}
		};
		inline bool getPaintSpans()
			{	return paintSpans;	};

		inline void setPaintRealTime()
		{
			if(paintRealTime)
			{	paintRealTime = false;	}
			else
			{	paintRealTime = true;	}
		};

		inline bool isPaintRealTime()
		{
			return paintRealTime;
		}

		inline void setNumUsagesToDisplay(int n)
			{	numUsagesToDisplay = n;	};
		inline int getNumUsagesToDisplay()
			{	return numUsagesToDisplay;	};
		inline void setUsageHistStartTime(int t)
			{	usageHistStartTime = t;	};
		inline int getUsageHistStartTime()
			{	return usageHistStartTime;	};
		
		bool terminateProgram;

		inline void setRouting(const char* rout)
		{	sprintf(routing,"%s",rout);	};
		inline void setTopology(const char* topo)
		{	sprintf(topology,"%s",topo);	};
		inline void setWavelength(const char* wave)
		{	sprintf(wavelength,"%s",wave);	};
		inline void setProbing(const char* prob)
		{	sprintf(probing,"%s",prob);	};
		inline void setQualityAware(int q)
		{	qualityAware = (bool)q;	};

		inline void setOverallBlocking(double b)
		{	overallBlocking = b;	};
		inline void setCollisions(double c)
		{	collisions = c;	};
		inline void setBadQuality(double b)
		{	badquality = b;	};
		inline void setNonResource(double n)
		{	nonresource = n;	};
		inline void setAvgProbesPerRequest(double a)
		{	avgprobesperrequest = a;	};
		inline void setAvgRequestDelayTime(double t)
		{	avgrequestdelaytime = t;	};
		inline void setAvgConnHopCount(double a)
		{	avgconnhopcount = a;	};
		inline void setAvgConnSpanCount(double a)
		{	avgconnspancount = a;	};
		inline void setAvgASEnoise(double a)
		{	avgASEnoise = a;	}
		inline void setAvgFWMnoise(double a)
		{	avgFWMnoise = a;	};
		inline void setAvgXPMnoise(double a)
		{	avgXPMnoise = a;	};
		
		void detailScreen();
#endif

		inline size_t getNumberOfConnections()
			{ return numberOfConnections; };

		inline double generateRandomZeroToOne()
			{ return generateZeroToOne(generator); }

		std::default_random_engine generator;

		std::uniform_real_distribution<double> generateZeroToOne;
		std::uniform_int_distribution<size_t> generateRandomRouter;

		std::exponential_distribution<double> generateRandomDuration;
		std::exponential_distribution<double> generateArrivalInterval;

		inline const std::string& getTopology()
			{ return topology; };

		inline void setTopology(const std::string& t)
			{ topology = t; };

		inline MessageLogger* getLogger()
			{ return logger; };

		static const int MORE_SIMULATIONS;
		static const int COMPLETED_ALL_SIMULATIONS;

	private:
		std::vector<Router*> routers;
		std::vector<Workstation*> workstations;

		void setTopologyParameters(const std::string& f);
		void setWorkstationParameters(const std::string& f);

		EventQueue* queue;

		double globalTime;

		MessageLogger* logger;

		ResourceManager* rm;

		std::string RoutingAlgorithmNames[NUMBER_OF_ROUTING_ALGORITHMS];
		std::string WavelengthAlgorithmNames[NUMBER_OF_WAVELENGTH_ALGORITHMS];
		std::string ProbeStyleNames[NUMBER_OF_PROBE_STYLES];

		void setAlgorithmParameters(const std::string& f, unsigned int iterationCount);

		void activate_workstations();
		void deactivate_workstations();
		void generateTrafficEvent(size_t session);

		void update_link_usage();

#ifdef RUN_GUI
		void update_gui();
#endif
		void connection_request(ConnectionRequestEvent* cre);
		void create_connection_probe(CreateConnectionProbeEvent* ccpe);
		void create_connection_confirmation(CreateConnectionConfirmationEvent* ccce);
		void destroy_connection_probe(DestroyConnectionProbeEvent* dcpe);
		void collision_notification(CollisionNotificationEvent* cne);

		GlobalStats stats;

		bool isLoadPrevious;

		size_t numOfWavelengths;

#ifdef RUN_GUI
		double maxMaxUsage;
		bool paintDir; //true if all edges must be repainted.
		bool paintSpans;
		bool paintRealTime;
		int progBarLengthpx; //length of full progress bar
		double ConnsPerPx;
		double multFactor;
		int paintedConns;
		int numUsagesToDisplay;
		int usageHistStartTime;
		char topoFile[100];
		char wkstFile[100];//TODO: do I ever use these?

		//SUMMARY variables to be saved to file
		char topology[20];//TODO: for this group of variables,
		char routing[20];//there may be other variables already 
		char wavelength[20];//fulfilled which serve the same purpose.
		char probing[20];//if possible, use those instead and
		bool qualityAware;//remember to get rid of all the set
		int totalActiveWKS;//methods used here
		int availableWaves;
		double overallBlocking;
		double collisions;
		double badquality;
		double nonresource;
		double avgprobesperrequest;
		double avgrequestdelaytime;
		double avgconnhopcount;
		double avgconnspancount;
		double avgASEnoise;
		double avgFWMnoise;
		double avgXPMnoise;

#endif

		unsigned int controllerIndex;

		RoutingAlgorithm CurrentRoutingAlgorithm;
		WavelengthAlgorithm CurrentWavelengthAlgorithm;
		ProbeStyle CurrentProbeStyle;
		bool CurrentQualityAware;
		unsigned int CurrentActiveWorkstations;

		QualityParameters qualityParams;
		void setQualityParameters(const std::string& f);

		size_t randomSeed;

		size_t numberOfRouters;
		size_t numberOfWorkstations;
		size_t numberOfEdges;

		size_t numberOfConnections;

		bool order_init;
		size_t* workstationOrder;

		double calculateDelay(size_t spans);

		bool sendResponse(CreateConnectionProbeEvent* probe);
		void clearResponses(CreateConnectionProbeEvent* probe);
		bool moreProbes(CreateConnectionProbeEvent* probe);
		int otherResponse(CreateConnectionProbeEvent* probe);

		void updateQMDegredation(Edge **connectionPath, size_t connectionLength, long long int wavelength);
		void updateQFactorStats(Edge **connectionPath, size_t connectionLength, long long int wavelength);

		double minDuration;

		size_t maxSpans;

		kShortestPathReturn* getKShortestPaths(ConnectionRequestEvent *cre, size_t probesToSend);
		CreateConnectionProbeEvent** calcProbesToSend(ConnectionRequestEvent *cre, kShortestPathReturn *kPath, 
			long long int &probesToSend, long long int &probeStart, long long int &probesSkipped);
		void sendProbes(ConnectionRequestEvent *cre, kShortestPathReturn *kPath, CreateConnectionProbeEvent** probesList,
			long long int probesToSend, long long int probeStart, long long int probesSkipped);

		std::vector<std::string> split(const std::string& s, char delimiter);

		int runCount;
		int maxRunCount;

		std::string topology;

		static const double TEN_HOURS;
		static const double SPEED_OF_LIGHT;
};

#endif
