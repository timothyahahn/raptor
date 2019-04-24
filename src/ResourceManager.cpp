// ____________________________________________________________________________
//
//  General Information:
//
//  File Name:      ResourceManager.cpp
//  Author:         Timothy Hahn, PhD
//  Project:        raptor
//
//  Description:    The file contains the implementation of the ResourceManager class.
//					The purpose of the ResourceManager is to calculate the path
//					from source to destination, calculate the Q-factor, and 
//					track the availability of the network resources.
//
//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//  Revision History:
//
//  05/20/2009	v1.0	Initial Version.
//  06/02/2009	v1.02	Minor optimizations and bug fixes.
//  04/14/2019  v2.0    Reworked version based upon cmake and octave
//
// ____________________________________________________________________________

#include "ErrorCodes.h"
#include "OctaveWrapper.h"
#include "ResourceManager.h"
#include "Thread.h"

#define HAVE_STRUCT_TIMESPEC
#include "pthread.h"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <sstream>

extern Thread* threadZero;
extern Thread** threads;

extern "C" void calc_k_shortest_paths(const kShortestPathParms &params, kShortestPathReturn* retVal);

///////////////////////////////////////////////////////////////////
//
// Function Name:	ResourceManager
// Description:		Default constructor with no arguments.
//
///////////////////////////////////////////////////////////////////
ResourceManager::ResourceManager() :
	fwm_combinations(nullptr), fwm_fs(nullptr), inter_indecies(nullptr), span_distance(nullptr), sys_fs_num(0),
	SP_paths(nullptr), kSP_edgeList(nullptr), wave_ordering(nullptr), sys_fs(new double[threadZero->getNumberOfWavelengths()]),
	sys_link_xpm_database(new double[threadZero->getNumberOfWavelengths() * threadZero->getNumberOfWavelengths()])
{
	calc_min_spans();

	OctaveWrapper::build_nonlinear_datastructure(sys_fs, sys_link_xpm_database);

	precompute_fwm_combinations();
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	~ResourceManager
// Description:		Default destructor with no arguments.
//
///////////////////////////////////////////////////////////////////
ResourceManager::~ResourceManager()
{
	delete[] sys_fs;
	delete[] sys_link_xpm_database;

	delete[] wave_ordering;

	delete[] span_distance;

	for(unsigned int s = 0; s < fwm_combinations->size(); ++s)
		fwm_combinations[s].clear();

	for(unsigned int w = 0; w < threadZero->getNumberOfWavelengths(); ++w)
	{
		delete[] (*fwm_fs)[w];
		delete[] (*inter_indecies)[w];
	}
	
	fwm_fs->clear();
	inter_indecies->clear();

	delete[] fwm_combinations;

	delete[] fwm_fs;
	delete[] inter_indecies;

	freeSPMatrix();

	delete[] kSP_edgeList;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	calculate_SP_path
// Description:		Calculates the shortest path from source to
//					destination for the SP algorithm
//
///////////////////////////////////////////////////////////////////
kShortestPathReturn* ResourceManager::calculate_SP_path(size_t src_index, size_t dest_index, size_t k, unsigned int ci)
{
	if(threads[ci]->getCurrentRoutingAlgorithm() == SHORTEST_PATH)
	{
		if(SP_paths != 0)
			if(SP_paths[src_index * threadZero->getNumberOfRouters() + dest_index] != 0)
				return SP_paths[src_index * threadZero->getNumberOfRouters() + dest_index];
	}

	if(kSP_edgeList == 0)
		build_KSP_EdgeList();

	kShortestPathParms kSP_params;

	kSP_params.src_node = src_index;
	kSP_params.dest_node = dest_index;
	kSP_params.k_paths = k;
	kSP_params.total_nodes = threadZero->getNumberOfRouters();
	kSP_params.total_edges = threadZero->getNumberOfEdges();
	kSP_params.edge_list = new kShortestPathEdges[kSP_params.total_edges];

	memcpy(kSP_params.edge_list,kSP_edgeList,sizeof(kShortestPathEdges) * kSP_params.total_edges);

	unsigned int num = 0;

	for(unsigned int a = 0; a < kSP_params.total_edges; ++a)
	{
		kSP_params.edge_list[a].edge_cost = 1;
	}

	kShortestPathReturn *kSP_return = new kShortestPathReturn();

	kSP_return->pathinfo = new size_t[kSP_params.k_paths * (kSP_params.total_nodes - 1)];
	kSP_return->pathcost = new double[kSP_params.k_paths];
	kSP_return->pathlen = new size_t[kSP_params.k_paths];

	calc_k_shortest_paths(kSP_params, kSP_return);

	delete[] kSP_params.edge_list;

	if(threads[ci]->getCurrentRoutingAlgorithm() == SHORTEST_PATH)
	{
		if(SP_paths != 0)
			SP_paths[src_index * threadZero->getNumberOfRouters() + dest_index] = kSP_return;
	}

	return kSP_return;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	calculate_span_distance
// Description:		Calculates the distance of the shortest path from 
//					source to destination for the SP algorithm.
//
///////////////////////////////////////////////////////////////////
unsigned int ResourceManager::calculate_span_distance(size_t src_index, size_t dest_index)
{
	unsigned int retVal = 0;

	if(kSP_edgeList == 0)
		build_KSP_EdgeList();

	kShortestPathParms kSP_params;

	kSP_params.src_node = src_index;
	kSP_params.dest_node = dest_index;
	kSP_params.k_paths = 1;
	kSP_params.total_nodes = threadZero->getNumberOfRouters();
	kSP_params.total_edges = threadZero->getNumberOfEdges();
	kSP_params.edge_list = new kShortestPathEdges[kSP_params.total_edges];

	memcpy(kSP_params.edge_list,kSP_edgeList,sizeof(kShortestPathEdges) * kSP_params.total_edges);

	unsigned int num = 0;

	kShortestPathReturn *kSP_return = new kShortestPathReturn();

	kSP_return->pathinfo = new size_t[1 * (kSP_params.total_nodes - 1)];
	kSP_return->pathcost = new double[1];
	kSP_return->pathlen = new size_t[1];

	calc_k_shortest_paths(kSP_params, kSP_return);

	delete[] kSP_params.edge_list;

	retVal = static_cast<unsigned int>(kSP_return->pathcost[0]);

	delete[] kSP_return->pathinfo;
	delete[] kSP_return->pathcost;
	delete[] kSP_return->pathlen;

	delete kSP_return;

	return retVal;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	calculate_LORA_path
// Description:		Calculates the shortest path from source to
//					destination for the LORA algorithm
//
///////////////////////////////////////////////////////////////////
kShortestPathReturn* ResourceManager::calculate_LORA_path(size_t src_index, size_t dest_index, size_t k, unsigned int ci)
{
	if(kSP_edgeList == 0)
		build_KSP_EdgeList();

	kShortestPathParms kSP_params;

	kSP_params.src_node = src_index;
	kSP_params.dest_node = dest_index;
	kSP_params.k_paths = k;
	kSP_params.total_nodes = threadZero->getNumberOfRouters();
	kSP_params.total_edges = threadZero->getNumberOfEdges();
	kSP_params.edge_list = new kShortestPathEdges[kSP_params.total_edges];

	memcpy(kSP_params.edge_list,kSP_edgeList,sizeof(kShortestPathEdges) * kSP_params.total_edges);

	unsigned int num = 0;

	for(unsigned int a = 0; a < threadZero->getNumberOfRouters(); ++a)
	{
		Router* routerA = threads[ci]->getRouterAt(a);

		for(unsigned int b = 0; b < threadZero->getNumberOfRouters(); ++b)
		{
			int edgeID = routerA->isAdjacentTo(b);

			if(edgeID >= 0)
			{
				kSP_params.edge_list[num].edge_cost = pow(threadZero->getBeta(),
					double(routerA->getEdgeByIndex(edgeID)->getAlgorithmUsage()));

				++num;
			}
		}
	}

	kShortestPathReturn *kSP_return = new kShortestPathReturn();

	kSP_return->pathinfo = new size_t[kSP_params.k_paths * (kSP_params.total_nodes - 1)];
	kSP_return->pathcost = new double[kSP_params.k_paths];
	kSP_return->pathlen = new size_t[kSP_params.k_paths];

	calc_k_shortest_paths(kSP_params, kSP_return);

	delete[] kSP_params.edge_list;

	return kSP_return;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	calculate_IA_path
// Description:		Calculates the shortest path from source to
//					destination for the IA-BF algorithm
//
///////////////////////////////////////////////////////////////////
kShortestPathReturn* ResourceManager::calculate_IA_path(size_t src_index, size_t dest_index, unsigned int ci)
{
	kShortestPathParms* kSP_params = new kShortestPathParms[threadZero->getNumberOfWavelengths()];

	for(unsigned int w = 0; w < threadZero->getNumberOfWavelengths(); ++w)
	{
		kSP_params[w].src_node = src_index;
		kSP_params[w].dest_node = dest_index;
		kSP_params[w].k_paths = 1;
		kSP_params[w].total_nodes = threadZero->getNumberOfRouters();
		kSP_params[w].total_edges = 0;
	
		kSP_params[w].edge_list = new kShortestPathEdges[threadZero->getNumberOfEdges()];
	}

	for(unsigned int a = 0; a < threadZero->getNumberOfRouters(); ++a)
	{
		Router* routerA = threads[ci]->getRouterAt(a);

		for(unsigned int b = 0; b < threadZero->getNumberOfRouters(); ++b)
		{
			int edgeID = routerA->isAdjacentTo(b);

			if(edgeID >= 0)
			{
				Edge* edge = routerA->getEdgeByDestination(b);

				for(unsigned int w = 0; w < threadZero->getNumberOfWavelengths(); ++w)
				{
					if(edge->getStatus(w) == EDGE_FREE)
					{
						kSP_params[w].edge_list[kSP_params[w].total_edges].src_node = a;
						kSP_params[w].edge_list[kSP_params[w].total_edges].dest_node = b;

						kSP_params[w].edge_list[kSP_params[w].total_edges].edge_cost = 
							double(edge->getNumberOfSpans());

						++kSP_params[w].total_edges;
					}
				}
			}
		}
	}

	kShortestPathReturn *kSP_temp = new kShortestPathReturn;

	kSP_temp->pathinfo = new size_t[kSP_params[0].total_nodes - 1];
	kSP_temp->pathcost = new double[1];
	kSP_temp->pathlen = new size_t[1];

	kShortestPathReturn *kSP_return = new kShortestPathReturn;

	kSP_return->pathinfo = new size_t[(kSP_params[0].total_nodes - 1) * threadZero->getNumberOfWavelengths()];
	kSP_return->pathcost = new double[threadZero->getNumberOfWavelengths()];
	kSP_return->pathlen = new size_t[threadZero->getNumberOfWavelengths()];

	for(unsigned int w = 0; w < threadZero->getNumberOfWavelengths(); ++w)
	{
		calc_k_shortest_paths(kSP_params[w],kSP_temp);

		kSP_return->pathcost[w] = kSP_temp->pathcost[0];
		kSP_return->pathlen[w] = kSP_temp->pathlen[0];

		if(kSP_return->pathlen[w] != std::numeric_limits<int>::infinity())
			for(unsigned int p = 0; p < kSP_return->pathlen[w]; ++p)
				kSP_return->pathinfo[w * (kSP_params[0].total_nodes - 1) + p] = kSP_temp->pathinfo[p];

		delete[] kSP_params[w].edge_list;
	}

	delete[] kSP_temp->pathinfo;
	delete[] kSP_temp->pathcost;
	delete[] kSP_temp->pathlen;

	delete kSP_temp;

	delete[] kSP_params;

	return kSP_return;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	calculate_PAR_path
// Description:		Calculates the shortest path from source to
//					destination for the PAR algorithm
//
///////////////////////////////////////////////////////////////////
kShortestPathReturn* ResourceManager::calculate_PAR_path(size_t src_index, size_t dest_index, size_t k, unsigned int ci)
{
	unsigned int iterationCount = 1;
	unsigned int kPathsFound = 0;
	unsigned int *kPathsStatus;

	const unsigned int PATH_VALID = 1;
	const unsigned int PATH_TOO_LONG = 2;
	const unsigned int PATH_INFINITY = 3;

	const unsigned int MAX_ITERATIONS = 4;

	while(kPathsFound < k && iterationCount <= MAX_ITERATIONS)
	{
		unsigned int kProduct	= 0;

		if(k == 1)
		{
			kProduct = static_cast<unsigned int>(pow(double(k+1),iterationCount-1));
		}
		else
		{
			kProduct = static_cast<unsigned int>(pow(double(k),iterationCount));
		}

		kShortestPathReturn* lora_ksp = calculate_LORA_path(src_index,dest_index,kProduct,ci);

		kPathsFound = 0;
		kPathsStatus = new unsigned int[kProduct];

		for(unsigned int a = 0; a < kProduct; ++a)
		{
			if(lora_ksp->pathcost[a] == std::numeric_limits<double>::infinity())
			{
				kPathsStatus[a] = PATH_INFINITY;
			}
			else
			{
				size_t pathSpans = 0;

				for(size_t p = 0; p < lora_ksp->pathlen[a] - 1; ++p)
				{
					pathSpans += threadZero->getRouterAt(lora_ksp->pathinfo[a * (threadZero->getNumberOfRouters() - 1) + p])
						->getEdgeByDestination(lora_ksp->pathinfo[a * (threadZero->getNumberOfRouters() - 1) + p + 1])->getNumberOfSpans();
				}

				if(pathSpans > threadZero->getMaxSpans())
				{
					kPathsStatus[a] = PATH_TOO_LONG;
				}
				else
				{
					kPathsStatus[a] = PATH_VALID;
					++kPathsFound;

					if(kPathsFound == k)
						break;
				}
			}
		}

		//Success. We have found the k shortest paths that satisfy the constraints,
		//so we return them.
		if(kPathsFound == k)
		{
			kShortestPathReturn *kSP_return = new kShortestPathReturn();

			kSP_return->pathinfo = new size_t[k * (threadZero->getNumberOfRouters() - 1)];
			kSP_return->pathcost = new double[k];
			kSP_return->pathlen = new size_t[k];
	
			unsigned int kIndex = 0;

			for(unsigned int a = 0; a < k; ++a)
			{
				while(kPathsStatus[kIndex] != PATH_VALID)
					++kIndex;

				kSP_return->pathcost[a] = lora_ksp->pathcost[kIndex];
				kSP_return->pathlen[a] = lora_ksp->pathlen[kIndex];

				for(unsigned int b = 0; b < kSP_return->pathlen[a]; ++b)
				{
					kSP_return->pathinfo[a * (threadZero->getNumberOfRouters() - 1) + b] =
						lora_ksp->pathinfo[kIndex * (threadZero->getNumberOfRouters() - 1) + b];
				}

				++kIndex;
			}

			delete[] lora_ksp->pathcost;
			delete[] lora_ksp->pathinfo;
			delete[] lora_ksp->pathlen;

			delete lora_ksp;

			delete[] kPathsStatus;

			return kSP_return;
		}

		//If we haven't found k paths that satisify the constraints after MAX_ITERATIONS,
		//then we just return the c paths that do satisfy the constraints and the first
		//k - c paths that do not satisfy the constraints.
		if(iterationCount == MAX_ITERATIONS)
		{
			kShortestPathReturn *kSP_return = new kShortestPathReturn();

			kSP_return->pathinfo = new size_t[k * (threadZero->getNumberOfRouters() - 1)];
			kSP_return->pathcost = new double[k];
			kSP_return->pathlen = new size_t[k];
	
			unsigned int kIndex = 0;
				
			//Return the c paths that satisify the constraints
			for(unsigned int a = 0; a < kPathsFound; ++a)
			{
				while(kPathsStatus[kIndex] != PATH_VALID)
					++kIndex;

				kSP_return->pathcost[a] = lora_ksp->pathcost[kIndex];
				kSP_return->pathlen[a] = lora_ksp->pathlen[kIndex];

				for(unsigned int b = 0; b < kSP_return->pathlen[a]; ++b)
				{
					kSP_return->pathinfo[a * (threadZero->getNumberOfRouters() - 1) + b] =
						lora_ksp->pathinfo[kIndex * (threadZero->getNumberOfRouters() - 1) + b];
				}

				++kIndex;
			}

			kIndex = 0;

			//Return the first k - c paths that do not satisfy the constraints
			for(unsigned int c = kPathsFound; c < k; ++c)
			{
				while(kPathsStatus[kIndex] != PATH_TOO_LONG)
					++kIndex;

				kSP_return->pathcost[c] = lora_ksp->pathcost[kIndex];
				kSP_return->pathlen[c] = lora_ksp->pathlen[kIndex];

				for(unsigned int d = 0; d < kSP_return->pathlen[c]; ++d)
				{
					kSP_return->pathinfo[c * (threadZero->getNumberOfRouters() - 1) + d] =
						lora_ksp->pathinfo[kIndex * (threadZero->getNumberOfRouters() - 1) + d];
				}

				++kIndex;
			}

			delete[] lora_ksp->pathcost;
			delete[] lora_ksp->pathinfo;
			delete[] lora_ksp->pathlen;

			delete lora_ksp;

			delete[] kPathsStatus;

			return kSP_return;
		}

		++iterationCount;

		delete[] lora_ksp->pathcost;
		delete[] lora_ksp->pathinfo;
		delete[] lora_ksp->pathlen;

		delete lora_ksp;

		delete[] kPathsStatus;
	}

	return 0;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	calculate_QM_path
// Description:		Calculates the shortest path from source to
//					destination for the QM algorithm
//
///////////////////////////////////////////////////////////////////
kShortestPathReturn* ResourceManager::calculate_QM_path(size_t src_index, size_t dest_index, size_t k, unsigned int ci)
{
	if(kSP_edgeList == 0)
		build_KSP_EdgeList();

	kShortestPathParms kSP_params;

	kSP_params.src_node = src_index;
	kSP_params.dest_node = dest_index;
	kSP_params.k_paths = k;
	kSP_params.total_nodes = threadZero->getNumberOfRouters();
	kSP_params.total_edges = threadZero->getNumberOfEdges();
	kSP_params.edge_list = new kShortestPathEdges[kSP_params.total_edges];

	memcpy(kSP_params.edge_list,kSP_edgeList,sizeof(kShortestPathEdges) * kSP_params.total_edges);

	unsigned int num = 0;

	for(unsigned int a = 0; a < threadZero->getNumberOfRouters(); ++a)
	{
		Router* routerA = threads[ci]->getRouterAt(a);

		for(unsigned int b = 0; b < threadZero->getNumberOfRouters(); ++b)
		{
			int edgeID = routerA->isAdjacentTo(b);

			if(edgeID >= 0)
			{
				kSP_params.edge_list[num].edge_cost = routerA->getEdgeByIndex(edgeID)->getQMDegredation();

				++num;
			}
		}
	}

	kShortestPathReturn *kSP_return = new kShortestPathReturn();

	kSP_return->pathinfo = new size_t[kSP_params.k_paths * (kSP_params.total_nodes - 1)];
	kSP_return->pathcost = new double[kSP_params.k_paths];
	kSP_return->pathlen = new size_t[kSP_params.k_paths];

	calc_k_shortest_paths(kSP_params, kSP_return);

	delete[] kSP_params.edge_list;

	return kSP_return;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	calculate_AQoS_path
// Description:		Calculates the shortest path from source to
//					destination for the AQoS algorithm
//
///////////////////////////////////////////////////////////////////
kShortestPathReturn* ResourceManager::calculate_AQoS_path(size_t src_index, size_t dest_index, size_t k, unsigned int ci)
{
	//First calculate the 2k shortest paths via the QM method
	kShortestPathReturn* QM_paths = calculate_QM_path(src_index, dest_index,k * 2,ci);
	size_t* QM_paths_availability = new size_t[k * 2];

	//Determine the amount of wavelengths available on each edge.
	for(unsigned int p = 0; p < k * 2; ++p)
	{
		if(QM_paths->pathcost[p] == std::numeric_limits<double>::infinity())
		{
			QM_paths_availability[p] = 0;
		}
		else
		{
			QM_paths_availability[p] = threadZero->getNumberOfWavelengths();

			for(size_t w = 0; w < threadZero->getNumberOfWavelengths(); ++w)
			{
				for(size_t r = 0; r < QM_paths->pathlen[p] - 1; ++r)
				{
					size_t srcIndex = QM_paths->pathinfo[p * (threadZero->getNumberOfRouters() - 1) + r];
					size_t destIndex = QM_paths->pathinfo[p * (threadZero->getNumberOfRouters() - 1) + r + 1];

					if(threads[ci]->getRouterAt(srcIndex)->getEdgeByDestination(destIndex)->getStatus(w) == EDGE_USED)
					{
						--QM_paths_availability[p];
						break;
					}

				}
			}
		}
	}

	//Copy the paths into the kpaths structure, ordered via the number of available wavelengths
	kShortestPathReturn *kSP_return = new kShortestPathReturn();

	kSP_return->pathinfo = new size_t[k * (threadZero->getNumberOfRouters() - 1)];
	kSP_return->pathcost = new double[k];
	kSP_return->pathlen = new size_t[k];

	for(size_t a = 0; a < k; ++a)
	{
		size_t maxAvailable = 0;
		size_t maxIndex = k * 2;

		for(size_t b = 0; b < k * 2; ++b)
		{
			if(maxAvailable < QM_paths_availability[b])
			{
				maxAvailable = QM_paths_availability[b];
				maxIndex = b;
			}
		}

		if(maxIndex != k * 2)
		{
			kSP_return->pathcost[a] = QM_paths->pathcost[maxIndex];
			kSP_return->pathlen[a] = QM_paths->pathlen[maxIndex];

			for(unsigned int c = 0; c < threadZero->getNumberOfRouters() - 1; ++c)
			{
				kSP_return->pathinfo[a * (threadZero->getNumberOfRouters() - 1) + c] =
					QM_paths->pathinfo[maxIndex * (threadZero->getNumberOfRouters() - 1) + c];
			}

			QM_paths_availability[maxIndex] = 0;
		}
		else
		{
			kSP_return->pathcost[a] = std::numeric_limits<double>::infinity();
			kSP_return->pathlen[a] = std::numeric_limits<int>::infinity();

			for(unsigned int c = 0; c < threadZero->getNumberOfRouters() - 1; ++c)
			{
				kSP_return->pathinfo[a * (threadZero->getNumberOfRouters() - 1) + c] = std::numeric_limits<int>::infinity();
			}
		}
	}

	//Delete the memory created
	delete[] QM_paths_availability;
	
	delete[] QM_paths->pathcost;
	delete[] QM_paths->pathinfo;
	delete[] QM_paths->pathlen;

	delete QM_paths;

	return kSP_return;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	calculate_ACO_path
// Description:		Calculates the shortest path from source to
//					destination for the ACO algorithm
//
///////////////////////////////////////////////////////////////////
kShortestPathReturn* ResourceManager::calculate_ACO_path(size_t src_index, size_t dest_index, size_t k, unsigned int ci)
{
	double alpha = threadZero->getQualityParams().DP_alpha;

	double l_exp = (double(threadZero->getMaxSpans()) + double(span_distance[src_index * threadZero->getNumberOfRouters() + dest_index])) / 2.0;
	double Q_exp = 10.0 * log10(threadZero->getQualityParams().channel_power/sqrt(l_exp * threadZero->getQualityParams().ASE_perEDFA[threadZero->getQualityParams().halfwavelength]));

	kShortestPathReturn* kSP_return = new kShortestPathReturn();

	kSP_return->pathcost = new double[k];
	kSP_return->pathlen = new size_t[k];
	kSP_return->pathinfo = new size_t[k * threadZero->getNumberOfRouters() - 1];

	for(unsigned int k1 = 0; k1 < k; ++k1)
	{
		kSP_return->pathcost[k1] = 0.0;
		kSP_return->pathlen[k1] = 0;
	}

	for(unsigned n = 0; n < threadZero->getNumberOfRouters(); ++n)
	{
		for(unsigned int e = 0; e < threads[ci]->getRouterAt(n)->getNumberOfEdges(); ++e)
		{
			threads[ci]->getRouterAt(n)->getEdgeByIndex(e)->resetPheremone(ci,
				span_distance[src_index * threadZero->getNumberOfRouters() + dest_index]);
		}

		threads[ci]->getRouterAt(n)->generateACOProbabilities(dest_index);
	}

	for(unsigned int i = 0; i < threadZero->getQualityParams().MM_ACO_N_iter; ++i)
	{
		Ant *ants = new Ant[threadZero->getQualityParams().ACO_ants];

		for(unsigned int a = 0; a < threadZero->getQualityParams().ACO_ants; ++a)
		{
			ants[a].location = threads[ci]->getRouterAt(src_index);
			ants[a].pathlen = 0;
			ants[a].path = new Edge*[threadZero->getNumberOfRouters() - 1];

			while(ants[a].location != threads[ci]->getRouterAt(dest_index) && ants[a].pathlen + 1 < (int)threadZero->getNumberOfRouters() - 1)
			{
				Edge* new_edge = ants[a].location->chooseEdge(threads[ci]->generateRandomZeroToOne());

				ants[a].location = threads[ci]->getRouterAt(new_edge->getDestinationIndex());
				ants[a].path[ants[a].pathlen] = new_edge;
				++ants[a].pathlen;

				for(int p = 0; p < ants[a].pathlen - 1; ++p)
				{
					if(ants[a].path[ants[a].pathlen - 1]->getDestinationIndex() == ants[a].path[p]->getSourceIndex())
					{
						ants[a].location = threads[ci]->getRouterAt(src_index);
						ants[a].pathlen = 0;

						break;
					}
				}
			}

			bool *free = new bool[threadZero->getNumberOfWavelengths()];

			for(unsigned int w1 = 0; w1 < threadZero->getNumberOfWavelengths(); ++w1)
			{
				free[w1] = true;
			}

			size_t spans = 0;

			for(int e = 0; e < ants[a].pathlen; ++e)
			{
				for(unsigned int w2 = 0; w2 < threadZero->getNumberOfWavelengths(); ++w2)
				{
					if(ants[a].path[e]->getStatus(w2) != EDGE_FREE)
						free[w2] = false;
				}

				spans += ants[a].path[e]->getNumberOfSpans();
			}

			double bestQ = 0.0;
			double pathWeight = 0.0;

			for(unsigned int w3 = 0; w3 < threadZero->getNumberOfWavelengths(); ++w3)
			{
				double ase;
				double fwm;
				double xpm;
				double Q;

				if(free[w3] == true)
				{
					Q = threadZero->getResourceManager()->estimate_Q(w3,ants[a].path,ants[a].pathlen,&xpm,&fwm,&ase,ci);

					if(Q > bestQ)
					{
						bestQ = Q;
					}
				}
			}

			delete[] free;

			pathWeight = (1.0 - alpha) * (bestQ / Q_exp) + alpha * l_exp / double(spans);

			if(pathWeight > kSP_return->pathcost[k-1] && bestQ  >= threadZero->getQualityParams().TH_Q)
			{
				//Check for duplicates....we need to keep k distinct paths!
				bool uniqueK = false;

				for(unsigned int k0 = 0; k0 < k; ++k0)
				{
					uniqueK = false;

					for(int r = 0; r < ants[a].pathlen; ++r)
					{
						if(ants[a].path[r]->getSourceIndex() !=
							kSP_return->pathinfo[k0 * (threadZero->getNumberOfRouters() - 1) + r])
						{
							uniqueK = true;
							break;
						}
					}

					if(uniqueK == false)
					{
						break;
					}
				}

				//We need to add unique paths only.
				if(uniqueK == true)
				{
					if(threads[ci]->getCurrentRoutingAlgorithm() == MAX_MIN_ACO)
					{
						i = 0;
					}

					//Calculate where to insert into the dest_node structure
					size_t k1 = k - 1;
					size_t k2 = 0;

					while(k1 > 0 && pathWeight > kSP_return->pathcost[k1-1])
					{
						--k1;
					}

					k2 = k - 1;

					while(k2 > k1)
					{
						for(unsigned int n = 0; n < threadZero->getNumberOfRouters() - 1; ++n)
						{
							kSP_return->pathinfo[k2 * (threadZero->getNumberOfRouters() - 1) + n] =
								kSP_return->pathinfo[(k2-1) * (threadZero->getNumberOfRouters() - 1) + n];
						}

						kSP_return->pathlen[k2] = kSP_return->pathlen[k2-1];
						kSP_return->pathcost[k2] = kSP_return->pathcost[k2-1];

						--k2;
					}

					//Insert where appropriate
					for(int n = 0; n < ants[a].pathlen; ++n)
					{
						kSP_return->pathinfo[k1 * (threadZero->getNumberOfRouters() - 1) + n] = 
							ants[a].path[n]->getSourceIndex();
					}

					kSP_return->pathinfo[k1 * (threadZero->getNumberOfRouters() - 1) + ants[a].pathlen] =
						ants[a].path[ants[a].pathlen - 1]->getDestinationIndex();

					kSP_return->pathlen[k1] = ants[a].pathlen + 1;
					kSP_return->pathcost[k1] = pathWeight;
				}
			}
		}

		for(unsigned n = 0; n < threadZero->getNumberOfRouters(); ++n)
		{
			for(unsigned int e = 0; e < threads[ci]->getRouterAt(n)->getNumberOfEdges(); ++e)
			{
				threads[ci]->getRouterAt(n)->getEdgeByIndex(e)->evaporatePheremone(ci);
			}
		}

		for(unsigned int a2 = 0; a2 < threadZero->getQualityParams().ACO_ants; ++a2)
		{
			if(threads[ci]->getCurrentRoutingAlgorithm() == ACO)
			{
				if(ants[a2].pathlen > 0)
				{
					size_t spans = 0;

					for(size_t n1 = 0; n1 < ants[a2].pathlen; ++n1)
					{
						spans += ants[a2].path[n1]->getNumberOfSpans();
					}
					for(int n2 = 0; n2 < ants[a2].pathlen; ++n2)
					{
						ants[a2].path[n2]->addPheremone(spans,ci);
					}				
				}
			}

			delete[] ants[a2].path;
		}

		if(threads[ci]->getCurrentRoutingAlgorithm() == MAX_MIN_ACO)
		{
			if(kSP_return->pathlen[0] > 0)
			{
				size_t spans = 0;

				for(unsigned int n1 = 0; n1 < kSP_return->pathlen[0] - 1; ++n1)
				{
					spans += threads[ci]->getRouterAt(kSP_return->pathinfo[n1])->
						getEdgeByDestination(kSP_return->pathinfo[n1+1])->getNumberOfSpans();
				}
				for(unsigned int n2 = 0; n2 < kSP_return->pathlen[0] - 1; ++n2)
				{
					threads[ci]->getRouterAt(kSP_return->pathinfo[n2])->getEdgeByDestination(kSP_return->pathinfo[n2+1])->
						addPheremone(spans,ci);
				}
			}
		}

		for(unsigned n2 = 0; n2 < threadZero->getNumberOfRouters(); ++n2)
		{
			threads[ci]->getRouterAt(n2)->generateACOProbabilities(dest_index);
		}

		delete[] ants;
	}

	for(unsigned int k1 = 0; k1 < k; ++k1)
	{
		if(kSP_return->pathlen[k1] == 0 && threads[ci]->getCurrentRoutingAlgorithm() == ACO)
		{
			kSP_return->pathcost[k1] = std::numeric_limits<double>::infinity();
			kSP_return->pathlen[k1] = std::numeric_limits<int>::infinity();
		}
	}

	return kSP_return;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	calculate_ACO_path
// Description:		Calculates the shortest path from source to
//					destination for the Max Min ACO algorithm
//
///////////////////////////////////////////////////////////////////
kShortestPathReturn* ResourceManager::calculate_MM_ACO_path(size_t src_index, size_t dest_index, size_t k, unsigned int ci)
{
	kShortestPathReturn* kSP_return = new kShortestPathReturn();
	kShortestPathReturn** mmACO_iters = new kShortestPathReturn*[threadZero->getQualityParams().MM_ACO_N_reset + 1];

	kSP_return->pathcost = new double[k];
	kSP_return->pathlen = new size_t[k];
	kSP_return->pathinfo = new size_t[k * threadZero->getNumberOfRouters() - 1];

	for(unsigned int k1 = 0; k1 < k; ++k1)
	{
		kSP_return->pathcost[k1] = 0.0;
		kSP_return->pathlen[k1] = 0;
	}

	for(unsigned int r = 0; r <= threadZero->getQualityParams().MM_ACO_N_reset; ++r)
	{
		mmACO_iters[r] = this->calculate_ACO_path(src_index, dest_index,k,ci);

		for(unsigned int k2 = 0; k2 < k; ++k2)
		{
			if(mmACO_iters[r]->pathcost[k2] > kSP_return->pathcost[k-1])
			{
				//Calculate where to insert into the dest_node structure
				size_t k3 = k - 1;
				size_t k4 = 0;

				while(k3 > 0 && mmACO_iters[r]->pathcost[k2] > kSP_return->pathcost[k3-1])
				{
					--k3;
				}

				k4 = k - 1;

				while(k4 > k3)
				{
					for(unsigned int n = 0; n < threadZero->getNumberOfRouters() - 1; ++n)
					{
						kSP_return->pathinfo[k4 * (threadZero->getNumberOfRouters() - 1) + n] =
							kSP_return->pathinfo[(k4-1) * (threadZero->getNumberOfRouters() - 1) + n];
					}

					kSP_return->pathlen[k4] = kSP_return->pathlen[k4-1];
					kSP_return->pathcost[k4] = kSP_return->pathcost[k4-1];

					--k4;
				}

				//Insert where appropriate
				for(unsigned int n = 0; n < mmACO_iters[r]->pathlen[k2]; ++n)
				{
					kSP_return->pathinfo[k3 * (threadZero->getNumberOfRouters() - 1) + n] = 
						mmACO_iters[r]->pathinfo[k2 * (threadZero->getNumberOfRouters() - 1) + n];
				}

				kSP_return->pathlen[k3] = mmACO_iters[r]->pathlen[k2];
				kSP_return->pathcost[k3] = mmACO_iters[r]->pathcost[k2];
			}
		}
	}

	for(unsigned int k5 = 0; k5 < k; ++k5)
	{
		if(kSP_return->pathlen[k5] == 0)
		{
			kSP_return->pathcost[k5] = std::numeric_limits<double>::infinity();
			kSP_return->pathlen[k5] = std::numeric_limits<int>::infinity();
		}
	}

	for(unsigned int r = 0; r <= threadZero->getQualityParams().MM_ACO_N_reset; ++r)
	{
		delete[] mmACO_iters[r]->pathcost;
		delete[] mmACO_iters[r]->pathinfo;
		delete[] mmACO_iters[r]->pathlen;

		delete mmACO_iters[r];
	}

	delete[] mmACO_iters;

	return kSP_return;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	calculate_DP_path
// Description:		Calculates the shortest path from source to
//					destination for the DP algorithm
//
///////////////////////////////////////////////////////////////////
kShortestPathReturn* ResourceManager::calculate_DP_path(size_t src_index, size_t dest_index, size_t k, unsigned int ci)
{
	double alpha = threadZero->getQualityParams().DP_alpha;

	size_t origK = k;
	
	if(threadZero->getTopology().compare("NSF") == 0 && k > 1 && k < 4)
	{
		k = 4;
	}
	else if(threadZero->getTopology().compare("Mesh") == 0 && k > 1 && k < 7)
	{
		k = 7;
	}
	else if(threadZero->getTopology().compare("Mesh8x8") == 0 && k > 1 && k < 7)
	{
		k = 7;
	}

	double l_exp = (double(threadZero->getMaxSpans()) + double(span_distance[src_index * threadZero->getNumberOfRouters() + dest_index])) / 2.0;
	double Q_exp = 10.0 * log10(threadZero->getQualityParams().channel_power/sqrt(l_exp * threadZero->getQualityParams().ASE_perEDFA[threadZero->getQualityParams().halfwavelength]));

	for(unsigned int r1 = 0; r1 < threadZero->getNumberOfRouters(); ++r1)
	{
		DP_node* node = new DP_node();

		node->paths = new Edge*[k * (threadZero->getNumberOfRouters() - 1)];
		node->waveAvailability = new bool[k * threadZero->getNumberOfWavelengths()];
		node->pathLength = new unsigned int[k];
		node->pathSpans = new size_t[k];
		node->optimalWave = new unsigned int[k];
		node->pathQuality = new double[k];
		node->pathWeight = new double[k];

		for(unsigned int t = 0; t < k * (threadZero->getNumberOfRouters() - 1); ++t)
		{
			node->paths[t] = 0;
		}

		for(unsigned int w = 0; w < k * threadZero->getNumberOfWavelengths(); ++w)
		{
			node->waveAvailability[w] = true;
		}

		for(unsigned int k1 = 0; k1 < k; ++k1)
		{
			node->pathLength[k1] = 0;
			node->pathSpans[k1] = 0;
			node->pathQuality[k1] = 0.0;
			node->pathWeight[k1] = 0.0;
		}

		threads[ci]->getRouterAt(r1)->dp_node = node;
	}

	std::queue<DP_item*> Q;

	for(unsigned int e = 0; e < threads[ci]->getRouterAt(src_index)->getNumberOfEdges(); ++e)
	{
		bool addEdge = false;
		DP_item* item = new DP_item;

		item->path = new Edge*[threadZero->getNumberOfRouters() - 1];
		item->waveAvailability = new bool[threadZero->getNumberOfWavelengths()];

		item->path[0] = threads[ci]->getRouterAt(src_index)->getEdgeByIndex(e);

		item->pathLength = 1;
		item->pathSpans = item->path[0]->getNumberOfSpans();

		for(unsigned int n = 1; n < threadZero->getNumberOfRouters() - 1; ++n)
			item->path[n] = 0;

		for(unsigned int w = 0; w < threadZero->getNumberOfWavelengths(); ++w)
		{
			item->waveAvailability[w] = item->path[0]->getStatus(w) == EDGE_FREE;

			if(item->waveAvailability[w] == true)
				addEdge = true;
		}
		
		if(addEdge == true)
		{
			Q.push(item);
		}
		else
		{
			delete[] item->path;
			delete[] item->waveAvailability;

			delete item;
		}
	}

	while(Q.size() > 0)
	{
		DP_item *current_item = Q.front();
		Q.pop();

		Edge* edge = current_item->path[current_item->pathLength-1];

		DP_node *dest_node = threads[ci]->getRouterAt(edge->getDestinationIndex())->dp_node;

		unsigned int additionalSpans = span_distance[edge->getDestinationIndex() * threadZero->getNumberOfRouters() + dest_index];

		if(current_item->pathSpans + additionalSpans < threadZero->getMaxSpans())
		{
			//Search for the wavelength with the best weight
			double bestQ = 0.0;
			double pathWeight = 0.0;
			double waveWeight = 0.0;
			unsigned int bestW = 0;

			for(unsigned int w = 0; w < threadZero->getNumberOfWavelengths(); ++w)
			{
				double ase = 0.0;
				double xpm = 0.0;
				double fwm = 0.0;

				double Q = 0.0;

				double bestCaseQ = 0.0;
				double bestCaseASE = 0.0;

				if(current_item->waveAvailability[w] == true)
				{
					Q = threadZero->getResourceManager()->estimate_Q(w,&current_item->path[0],current_item->pathLength,&xpm,&fwm,&ase,ci);

					bestCaseASE = additionalSpans * threadZero->getQualityParams().ASE_perEDFA[threadZero->getQualityParams().halfwavelength];
					bestCaseQ = 10.0 * log10(threadZero->getQualityParams().channel_power/sqrt(bestCaseASE + ase + xpm + fwm));

					waveWeight = (1.0 - alpha) * (Q / Q_exp) + 
						alpha * l_exp / double(current_item->pathSpans + additionalSpans);

					if(waveWeight > pathWeight && Q > threadZero->getQualityParams().TH_Q)
					{
						bestQ = Q;
						pathWeight = waveWeight;
						bestW = w;

						if(alpha == 1.0)
						{
							break;
						}
					}
				}
			}

			if(bestQ >= threadZero->getQualityParams().TH_Q && pathWeight > dest_node->pathWeight[k - 1])
			{
				//Check for duplicates....we need to keep k distinct paths!
				bool uniqueK = false;

				for(unsigned int k0 = 0; k0 < k; ++k0)
				{
					uniqueK = false;

					for(unsigned int r = 0; r < current_item->pathLength; ++r)
					{
						if(current_item->path[r] != dest_node->paths[k0 * (threadZero->getNumberOfRouters() - 1) + r])
						{
							uniqueK = true;
							break;
						}
					}

					if(uniqueK == false)
					{
						break;
					}
				}

				//We need to add unique paths only.
				if(uniqueK != false)
				{
					//Calculate where to insert into the dest_node structure
					size_t k1 = k - 1;
					size_t k2 = 0;

					while(k1 > 0 && (pathWeight > dest_node->pathWeight[k1 - 1] || dest_node->pathLength[k1 -1] == 0))
					{
						--k1;
					}

					k2 = k - 1;

					while(k2 > k1)
					{
						memcpy(&dest_node->paths[k2 * (threadZero->getNumberOfRouters() - 1)],
							&dest_node->paths[(k2-1) * (threadZero->getNumberOfRouters() - 1)],
							sizeof(Edge*) * (threadZero->getNumberOfRouters() - 1));

						memcpy(&dest_node->waveAvailability[k2 * threadZero->getNumberOfWavelengths()],
							&dest_node->waveAvailability[(k2-1) * threadZero->getNumberOfWavelengths()],
							sizeof(bool) * threadZero->getNumberOfWavelengths());

						dest_node->pathLength[k2] = dest_node->pathLength[k2-1];
						dest_node->pathQuality[k2] = dest_node->pathQuality[k2-1];
						dest_node->optimalWave[k2] = dest_node->optimalWave[k2-1];
						dest_node->pathSpans[k2] = dest_node->pathSpans[k2-1];
						dest_node->pathWeight[k2] = dest_node->pathWeight[k2-1];

						--k2;
					}

					//Insert where appropriate
					memcpy(&dest_node->paths[k1 * (threadZero->getNumberOfRouters() - 1)],
						&current_item->path[0],	sizeof(Edge*) * (threadZero->getNumberOfRouters() - 1));

					memcpy(&dest_node->waveAvailability[k1 * threadZero->getNumberOfWavelengths()],
						&current_item->waveAvailability[0], sizeof(bool) * threadZero->getNumberOfWavelengths());

					dest_node->pathLength[k1] = current_item->pathLength;
					dest_node->pathSpans[k1] = current_item->pathSpans;
					dest_node->pathQuality[k1] = bestQ;
					dest_node->optimalWave[k1] = bestW;
					dest_node->pathWeight[k1] = pathWeight;

					for(unsigned int e = 0; e < threads[ci]->getRouterAt(edge->getDestinationIndex())->getNumberOfEdges(); ++e)
					{
						Edge* tmp_edge = threads[ci]->getRouterAt(edge->getDestinationIndex())->getEdgeByIndex(e);

						//We don't want to put edges with a destinaton of the source on the Q.
						//This would result in a cycle.
						if(tmp_edge->getDestinationIndex() == src_index && tmp_edge->getSourceIndex() == dest_index)
							continue;

						bool cycle = false;

						for(unsigned int r = 0; r < current_item->pathLength; ++r)
						{
							if(current_item->path[r]->getSourceIndex() == tmp_edge->getDestinationIndex())
							{
								cycle = true;
								break;
							}
						}

						if(cycle == true)
							continue;

						bool addEdge = false;
						DP_item* item = new DP_item;

						item->path = new Edge*[threadZero->getNumberOfRouters() - 1];
						item->waveAvailability = new bool[threadZero->getNumberOfWavelengths()];

						item->pathLength = current_item->pathLength + 1;
						item->pathSpans = current_item->pathSpans + tmp_edge->getNumberOfSpans();

						if(item->pathSpans > threadZero->getMaxSpans() ||
							(alpha == 0 && item->pathSpans > threadZero->getRouterAt(tmp_edge->getDestinationIndex())->dp_node->pathSpans[k-1] &&
							threadZero->getRouterAt(tmp_edge->getDestinationIndex())->dp_node->pathSpans[k-1] != 0))
						{
							delete[] item->path;
							delete[] item->waveAvailability;

							delete item;

							continue;
						}

						memcpy(item->path,current_item->path,sizeof(Edge*) * (threadZero->getNumberOfRouters() - 1));

						item->path[item->pathLength-1] = tmp_edge;

						for(unsigned int w = 0; w < threadZero->getNumberOfWavelengths(); ++w)
						{
							item->waveAvailability[w] = current_item->waveAvailability[w] && (tmp_edge->getStatus(w) == EDGE_FREE);

							if(item->waveAvailability[w] == true)
								addEdge = true;
						}
						
						if(addEdge == true)
						{
							Q.push(item);
						}
						else
						{
							delete[] item->path;
							delete[] item->waveAvailability;

							delete item;
						}
					}
				}
			}
		}

		delete[] current_item->path;
		delete[] current_item->waveAvailability;

		delete current_item;
	}

	k = origK;

	//Populate return structure with data from DP_node
	kShortestPathReturn* kSP_return = new kShortestPathReturn();

	kSP_return->pathcost = new double[k];
	kSP_return->pathlen = new size_t[k];
	kSP_return->pathinfo = new size_t[k * threadZero->getNumberOfRouters() - 1];

	DP_node* final_dp_node  = threads[ci]->getRouterAt(dest_index)->dp_node;

	for(unsigned int k1 = 0; k1 < k; ++k1)
	{
		if(final_dp_node->pathLength[k1] > 0)
		{
			kSP_return->pathcost[k1] = final_dp_node->optimalWave[k1];
			kSP_return->pathlen[k1] = final_dp_node->pathLength[k1] + 1;

			for(unsigned int r = 0; r < kSP_return->pathlen[k1] - 1; ++r)
			{
				kSP_return->pathinfo[k1 * (threadZero->getNumberOfRouters() - 1) + r] = 
					final_dp_node->paths[k1 * (threadZero->getNumberOfRouters() - 1) + r]->getSourceIndex();
			}

			kSP_return->pathinfo[k1 * (threadZero->getNumberOfRouters() - 1) + kSP_return->pathlen[k1] - 1] = 
				final_dp_node->paths[k1 * (threadZero->getNumberOfRouters() - 1) + kSP_return->pathlen[k1] - 2]->getDestinationIndex();
		}
		else
		{
			kSP_return->pathcost[k1] = std::numeric_limits<double>::infinity();
			kSP_return->pathlen[k1] = std::numeric_limits<int>::infinity();
		}
	}

	for(unsigned int r1 = 0; r1 < threadZero->getNumberOfRouters(); ++r1)
	{
		DP_node *node = threads[ci]->getRouterAt(r1)->dp_node;

		delete[] node->paths;
		delete[] node->waveAvailability;
		delete[] node->pathLength;
		delete[] node->pathSpans;
		delete[] node->optimalWave;
		delete[] node->pathQuality;
		delete[] node->pathWeight;

		delete node;
	}

	return kSP_return;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	choose_wavelength
// Description:		Calculates the shortest path from source to
//					destination for the SP algorithm
//
///////////////////////////////////////////////////////////////////
int ResourceManager::choose_wavelength(CreateConnectionProbeEvent* ccpe, unsigned int ci)
{
	int retval;

	bool* wave_available = new bool[threadZero->getNumberOfWavelengths()];
	size_t numberAvailableWaves = threadZero->getNumberOfWavelengths();

	for(unsigned int w = 0; w < threadZero->getNumberOfWavelengths(); ++w)
		wave_available[w] = true;

	for(unsigned int r = 0; r < ccpe->connectionLength; ++r)
	{
		Edge* edge = ccpe->connectionPath[r];

		for(size_t k = 0; k < threadZero->getNumberOfWavelengths(); ++k)
		{
			if(edge->getStatus(k) != EDGE_FREE && wave_available[k] == true)
			{
				wave_available[k] = false;
				--numberAvailableWaves;
			}
		}

	}

	if(numberAvailableWaves == 0)
	{
		delete[] wave_available;
		return NO_PATH_FAILURE;
	}

	if(threads[ci]->getCurrentWavelengthAlgorithm() == FIRST_FIT)
	{
		retval = first_fit(ccpe,ci,wave_available);
		delete[] wave_available;
	}
	else if(threads[ci]->getCurrentWavelengthAlgorithm() == FIRST_FIT_ORDERED)
	{
		retval = first_fit_with_ordering(ccpe,ci,wave_available);
		delete[] wave_available;
	}
	else if(threads[ci]->getCurrentWavelengthAlgorithm() == RANDOM_FIT)
	{
		retval = random_fit(ccpe,ci,wave_available,numberAvailableWaves);
		delete[] wave_available;
	}
	else if(threads[ci]->getCurrentWavelengthAlgorithm() == MOST_USED)
	{
		retval = most_used(ccpe,ci,wave_available);
		delete[] wave_available;
	}
	else if(threads[ci]->getCurrentWavelengthAlgorithm() == QUAL_FIRST_FIT)
	{
		return quality_first_fit(ccpe,ci,wave_available,numberAvailableWaves);
	}
	else if(threads[ci]->getCurrentWavelengthAlgorithm() == QUAL_FIRST_FIT_ORDERED)
	{
		return quality_first_fit_with_ordering(ccpe,ci,wave_available,numberAvailableWaves);
	}
	else if(threads[ci]->getCurrentWavelengthAlgorithm() == QUAL_RANDOM_FIT)
	{
		return quality_random_fit(ccpe,ci,wave_available,numberAvailableWaves);
	}
	else if(threads[ci]->getCurrentWavelengthAlgorithm() == QUAL_MOST_USED)
	{
		return quality_most_used(ccpe,ci,wave_available,numberAvailableWaves);
	}
	else if(threads[ci]->getCurrentWavelengthAlgorithm() == LEAST_QUALITY)
	{
		return least_quality_fit(ccpe,ci,wave_available);
	}
	else if(threads[ci]->getCurrentWavelengthAlgorithm() == MOST_QUALITY)
	{
		return most_quality_fit(ccpe,ci,wave_available);
	}
	else
	{
		threadZero->recordEvent("ERROR: Invalid value for CurrentWavelengthAlgorithm.\n",true,ci);
		exit(ERROR_CHOOSE_WAVELENGTH_2);
	}

	double q_factor = 0.0;
	double xpm_noise = 0.0;
	double fwm_noise = 0.0;
	double ase_noise = 0.0;

	if(retval >= 0)
	{
		q_factor = threadZero->getResourceManager()->estimate_Q(
			retval,ccpe->connectionPath,ccpe->connectionLength,&xpm_noise,
			&fwm_noise,&ase_noise,ci);

		if(threads[ci]->getCurrentQualityAware() == true && 
		   q_factor < threadZero->getQualityParams().TH_Q)
		{
			retval = QUALITY_FAILURE;
			ccpe->wavelength = QUALITY_FAILURE;
		}
	}

	print_connection_info(ccpe,q_factor,ase_noise,fwm_noise,xpm_noise,ci);

	return retval;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	estimate_Q
// Description:		Estimates the Q-factor based upon XPM noise,
//					FWM noise, and ASE noise
//
///////////////////////////////////////////////////////////////////
double ResourceManager::estimate_Q(int lambda, Edge **Path, size_t pathLen, double *xpm, double *fwm, double *ase, unsigned int ci)
{
	double noise = 0.0;
	double Q = 0.0;

	if(lambda >= 0 && lambda < static_cast<int>(threadZero->getNumberOfWavelengths()))
	{
		*xpm = path_xpm_noise(lambda, Path, pathLen, ci);
		*fwm = path_fwm_noise(lambda, Path, pathLen, ci);
	}
	else
	{
		*xpm = 0.0;
		*fwm = 0.0;
	}

	*ase = path_ase_noise(lambda, Path, pathLen, ci);

	noise = *xpm + *fwm + *ase;
  
	Q = 10.0 * log10(threadZero->getQualityParams().channel_power/sqrt(noise));

	return Q;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	path_ase_noise
// Description:		Calculates the ASE noise created along the path.
//
///////////////////////////////////////////////////////////////////
double ResourceManager::path_ase_noise(int lambda, Edge **Path, size_t pathLen, unsigned int ci)
{
	double spans = 0.0;
	
	for(unsigned int r = 0; r < pathLen; ++r)
	{
		spans += Path[r]->getNumberOfSpans();
	}

	return spans * threadZero->getQualityParams().ASE_perEDFA[lambda];
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	path_xpm_noise
// Description:		Calculates the ASE noise created along the path.
//
///////////////////////////////////////////////////////////////////
double ResourceManager::path_xpm_noise(int lambda, Edge **Path, size_t pathLen, unsigned int ci)
{
 	double noise = 0.0;  

	for(int wave = 0; wave < static_cast<int>(threadZero->getNumberOfWavelengths()); ++wave)
	{
		//We don't want to compute the XPM for cases where the wavelength is outside of the halfwin window
		//or where the wavelength is equal to the connection wavelength.
		if ((abs(wave - lambda) > threadZero->getQualityParams().nonlinear_halfwin) || wave == lambda) 
			continue;

		unsigned int index = 0; 
		size_t path_len = 0;
		
		while(index < pathLen)
		{
			unsigned int j = 0;

			for(j = index; j < pathLen; ++j)
			{
				if(Path[j]->getStatus(wave) == EDGE_USED)
				{
					//If the cumulative path length is zero, then we don't care about the session.
					if(path_len == 0)
					{
						path_len = Path[j]->getNumberOfSpans();
					}
					//If the path_len is greater than zero, then we need to make sure that the session number is
					//equal to the previous session.
					else if(Path[j-1]->getActiveSession(wave) == Path[j]->getActiveSession(wave))
					{
						path_len += Path[j]->getNumberOfSpans();
					}
					//If the link is active but the sessions are different, we need to consider this link on the
					//next iteration. Thus we need to break from the for loop.
					else
					{
						break;
					}
				}
				else
				{
					break;
				}
			}

			if(path_len == 0)
				++index;
			else
				index = j;

			if (path_len > 0)
				noise += path_xpm_term(path_len, lambda, wave);
		 
			path_len = 0;
			
		}
	}

	return noise;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	path_xpm_term
// Description:		Loads the XPM term from the nonlinear database
//
///////////////////////////////////////////////////////////////////
double ResourceManager::path_xpm_term(size_t spans, size_t lambda, size_t wave)
{
	return sys_link_xpm_database[lambda * threadZero->getNumberOfWavelengths() + wave] * double(spans) * double(spans);
}

///////////////////////////////////////////////////////////////////
//
// Function Name:    path_fwm_noise
// Description:      Calculates the ASE noise created along the path.
//
///////////////////////////////////////////////////////////////////
double ResourceManager::path_fwm_noise(int lambda, Edge **Path, size_t pathLen, unsigned int ci)
{
    double noise = 0.0;

	for(int r = 0; r < static_cast<int>(fwm_combinations[lambda].size() / 4); r++)
    {
        int i_id = fwm_combinations[lambda][r * 4 + 0];
        int j_id = fwm_combinations[lambda][r * 4 + 1];
        int k_id = fwm_combinations[lambda][r * 4 + 2];
        int d = fwm_combinations[lambda][r * 4 + 3];

		int i_wave = (*inter_indecies)[lambda][i_id];
        int j_wave = (*inter_indecies)[lambda][j_id];
        int k_wave = (*inter_indecies)[lambda][k_id];
        double fi = (*fwm_fs)[lambda][i_id];
        double fj = (*fwm_fs)[lambda][j_id];
        double fk = (*fwm_fs)[lambda][k_id];       

        unsigned int index = 0;
        size_t plen = 0;
		unsigned int j = 0;

        while(index < pathLen)
	    {
		    for(j = index; j < pathLen; ++j)
			{
				if((Path[j]->getStatus(i_wave) == EDGE_USED || i_wave == lambda) &&
	               (Path[j]->getStatus(j_wave) == EDGE_USED || j_wave == lambda) &&
	               (Path[j]->getStatus(k_wave) == EDGE_USED || k_wave == lambda))
				{
					//If the cumulative path length is zero, then we don't care about the session.
					if(plen == 0)
					{
						plen = Path[j]->getNumberOfSpans();
					}
					//If the path_len is greater than zero, then we need to make sure that the session numbers are
					//equal to the previous session.
					else if(Path[j-1]->getActiveSession(i_wave) == Path[j]->getActiveSession(i_wave) &&
							Path[j-1]->getActiveSession(j_wave) == Path[j]->getActiveSession(j_wave) &&
							Path[j-1]->getActiveSession(k_wave) == Path[j]->getActiveSession(k_wave))
					{
						plen += Path[j]->getNumberOfSpans();
					}
					//If the link is active but the sessions are different, we need to consider this link on the
					//next iteration. Thus we need to break from the for loop.
					else
					{
						break;
					}
				}
				else
				{
					break;
				}
			}
			
			index = j;
        
			if(plen > 0)
			{
				noise += path_fwm_term(plen,fi,fj,fk,sys_fs[lambda],d);
			}
			else if(plen == 0)
			{
				++index;
			}
        
			plen = 0;
		}//end while     
    } //end for
  
    return 2.0 * threadZero->getQualityParams().channel_power * noise;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	build_FWM_fs
// Description:		Generates the frequency combinations and returns
//					the number.
//
///////////////////////////////////////////////////////////////////
int ResourceManager::build_FWM_fs(double *inter_fs,int *inter_indecies, int lambda)
{
	int num = 0;
 
	for(unsigned int i = 0; i < threadZero->getNumberOfWavelengths(); ++i)
		if (abs(static_cast<int>(i - lambda)) <= threadZero->getQualityParams().nonlinear_halfwin)
		{
        	inter_indecies[num] = i;
			inter_fs[num] = sys_fs[i];
			num++;
		}
	
	return num;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	wave_combines
// Description: Return all FWM terms which can generate fc. The fs_coms
//				include indicies of fi, fj, fk and d. fs should be sorted.
//				Incicies are relative to fs.
//
///////////////////////////////////////////////////////////////////
int ResourceManager::wave_combines(double fc, double *fs,int fs_num,std::vector<int> &fs_coms)
{
	int num = 0;

	for(int i = 0; i < fs_num; i++)
		for(int j=0; j < fs_num; j++)
			for(int k = 0; k < fs_num; k++)
			{
				int fi = -1;
				int fj = -1;
				int fk = -1;
				int d =-1;

				if (-fs[i]+fs[j]+fs[k] == fc)
				{
					fi = j;
					fj = k; 
					fk = i; 
					d = degeneracy(fi,fj,fk);
				}
				else if (fs[i]-fs[j]+fs[k] == fc)
				{
					fi = i;
					fj = k;
					fk = j; 
					d = degeneracy(fi,fj,fk);
			    }
				else if (fs[i]+fs[j]-fs[k] == fc)
				{
					fi = i;
					fj = j;
					fk = k; 
					d = degeneracy(fi,fj,fk);
			    }
				else
				{
					continue;
				}

				if(d == -1)
				{
					continue;
				}

				int tmp1 = fi;
				int tmp2 = fj;
				
				if(tmp1 > tmp2)
				{
					fi = tmp2;
					fj = tmp1;
				}
				else
				{
					fi = tmp1;
					fj = tmp2;
				}
			
				if (!can_find(fi,fj,fk,fs_coms,num))
				{
					fs_coms.push_back(fi);   // fi, fj are sorted
					fs_coms.push_back(fj);
					fs_coms.push_back(fk);
					fs_coms.push_back(d);
					++num;
				}
			}

	return num;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	can_find
// Description: Decide whether a combination exists in the database.
//
///////////////////////////////////////////////////////////////////
bool ResourceManager::can_find(int fi,int fj,int fk,std::vector<int> &fs_coms,int com_num)
{
	for (int i=0;i<com_num;i++) 
		if (fi == fs_coms[i * 4 + 0] && fj == fs_coms[i * 4 + 1] && fk == fs_coms[i * 4 + 2])
			return true;
       
	return false;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	can_find
// Description: Decide the degeneracy factor of the FWM terms.
//
///////////////////////////////////////////////////////////////////
int ResourceManager::degeneracy(int fi,int fj,int fk)
{
	if (fi == fk || fj == fk) 
		return -1 ;            // excluding spm and xpm terms
	else
		if (fi==fj)
			return 3;
		else
			return 6;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	path_fwm_term
// Description: 
//
///////////////////////////////////////////////////////////////////
double ResourceManager::path_fwm_term(size_t spans,double fi,double fj, double fk,double fc,int dgen)
{ 
	double c = 2.99792457778e+8;
	double pi = 3.14159265358979323846;
	double lambdac = c / fc;
	double noise = 0.0;   
	double wdm_factor = 0.0;
	double Pi_0, Pj_0, Pk_0;

	double alpha = threadZero->getQualityParams().alpha;
	double channel_power = threadZero->getQualityParams().channel_power;
	double D = threadZero->getQualityParams().D;
	double S = threadZero->getQualityParams().S;
	double L = threadZero->getQualityParams().L;
	double gamma = threadZero->getQualityParams().gamma;
       
	if (fi==fc)
		Pi_0 = channel_power;
	else 
		Pi_0 = 0.5 * channel_power;
   
	if (fj==fc)
		Pj_0 = channel_power;
	else 
		Pj_0 = 0.5 * channel_power;
   
	if (fk==fc)
		Pk_0 = channel_power;
	else 
		Pk_0 = 0.5 * channel_power;
   
	double diff_kappa = 2.0 * pi * lambdac * lambdac / c * (fi - fc) * (fj - fc) * (D - lambdac * lambdac / c * (fi / 2.0 + fj / 2.0 - fc ) * S); 
	double diff_phi = 2.0 * pi * lambdac * lambdac / c * (fi - fc) * (fj - fc) * (-lambdac * lambdac / c * (fi / 2.0 + fj / 2.0 - fc) * S) * L;
	double Leff_square = (1.0 + exp(-2.0 * alpha * L) -  2.0 * exp(-alpha * L) * cos(diff_kappa * L) ) / (alpha * alpha + diff_kappa * diff_kappa);
	double tmp = gamma * gamma * dgen * dgen / 9.0 * Pi_0 * Pj_0 * Pk_0 * Leff_square;

    if (cos(diff_phi) != 1)
		wdm_factor = (double(1.0)-cos(diff_phi*spans)) /(double(1.0)-cos(diff_phi));
    else 
		wdm_factor  =  double(spans * spans);
           
    return tmp * wdm_factor;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	first_fit_wave
// Description:		Chooses a wavelength based upon the first fit
//					algorithm
//
///////////////////////////////////////////////////////////////////
int ResourceManager::first_fit(CreateConnectionProbeEvent* ccpe, unsigned int ci, bool* wave_available)
{
	for(unsigned int w = 0; w < threadZero->getNumberOfWavelengths(); ++w)
	{
		if(wave_available[w] == true)
		{
			return w;
		}
	}
	return NO_PATH_FAILURE;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	first_fit_with_ordering
// Description:		Chooses a wavelength based upon the first fit
//					with_ordering_algorithm
//
///////////////////////////////////////////////////////////////////
int ResourceManager::first_fit_with_ordering(CreateConnectionProbeEvent* ccpe, unsigned int ci, bool* wave_available)
{
	if(wave_ordering == 0)
	{
		generateWaveOrdering();
	}

	for(unsigned int w = 0; w < threadZero->getNumberOfWavelengths(); ++w)
	{
		if(wave_available[wave_ordering[w]] == true)
		{
			return wave_ordering[w];
		}
	}
	return NO_PATH_FAILURE;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	most_used
// Description:		Chooses a wavelength based upon the first fit
//					with_ordering_algorithm
//
///////////////////////////////////////////////////////////////////
int ResourceManager::most_used(CreateConnectionProbeEvent* ccpe, unsigned int ci, bool* wave_available)
{
	int *wave_counts = new int[threadZero->getNumberOfWavelengths()];

	for(unsigned int w = 0; w < threadZero->getNumberOfWavelengths(); ++w)
	{
		wave_counts[w] = 0;
	}

	for(unsigned int r = 0; r < threadZero->getNumberOfRouters(); ++r)
	{
		for(unsigned int e = 0; e < threads[ci]->getRouterAt(r)->getNumberOfEdges(); ++e)
		{
			Edge *edge = threads[ci]->getRouterAt(r)->getEdgeByIndex(e);

			for(unsigned int w = 0; w < threadZero->getNumberOfWavelengths(); ++w)
			{
				if(wave_available[w] == true)
				{
					if(edge->getStatus(w) == EDGE_USED)
					{
						++wave_counts[w];
					}
				}
			}
		}
	}

	int return_val = NO_PATH_FAILURE;
	int maxUsed = -1;

	for(unsigned int w = 0; w < threadZero->getNumberOfWavelengths(); ++w)
	{
		if(wave_counts[w] > maxUsed && wave_available[w] == true)
		{
			return_val = w;
			maxUsed = wave_counts[w];
		}
	}

	delete[] wave_counts;

	return return_val;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	random_fit
// Description:		Chooses a wavelength based upon the random
//					algorithm
//
///////////////////////////////////////////////////////////////////
int ResourceManager::random_fit(CreateConnectionProbeEvent* ccpe, unsigned int ci, bool* wave_available, size_t numberAvailableWaves)
{
	std::default_random_engine generator = std::default_random_engine(threadZero->getRandomSeed() * ccpe->sourceRouterIndex * ccpe->destinationRouterIndex * numberAvailableWaves);
	std::uniform_int_distribution<size_t> generateWavelength = std::uniform_int_distribution<size_t>(0, numberAvailableWaves - 1);

	size_t waveToReturn = generateWavelength(generator);
	size_t waveIndex = 0;

	for(size_t w = 0; w < threadZero->getNumberOfWavelengths(); ++w)
	{
		if(wave_available[w] == true)
		{
			if(waveIndex == waveToReturn)
			{
				return w;
			}
			else
			{
				++waveIndex;
			}
		}
	}

	threadZero->recordEvent("ERROR: Unexpected point in choose wavelength.\n",true,ci);
	exit(ERROR_CHOOSE_WAVELENGTH_1);
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	quality_first_fit
// Description:		Chooses a wavelength based upon the first fit
//					algorithm using quality
//
///////////////////////////////////////////////////////////////////
int ResourceManager::quality_first_fit(CreateConnectionProbeEvent* ccpe, unsigned int ci, bool* wave_available,size_t numberAvailableWaves)
{
	while(numberAvailableWaves > 0)
	{
		double xpm = 0.0;
		double fwm = 0.0;
		double ase = 0.0; 
		double Q_factor = 0.0;

		int wave = first_fit(ccpe,ci,wave_available);

		Q_factor = threadZero->getResourceManager()->estimate_Q(wave,ccpe->connectionPath,ccpe->connectionLength,&xpm,&fwm,&ase,ci);

		if(threadZero->getQualityParams().TH_Q <= Q_factor)
		{
			ccpe->wavelength = wave;

			print_connection_info(ccpe,Q_factor,ase,fwm,xpm,ci);

			delete[] wave_available;
			return wave;
		}
		else
		{
			--numberAvailableWaves;
			wave_available[wave] = false;
		}
	}

	delete[] wave_available;

	ccpe->wavelength = QUALITY_FAILURE;

	print_connection_info(ccpe,0.0,0.0,0.0,0.0,ci);

	return QUALITY_FAILURE;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	quality_first_fit_with_ordering
// Description:		Chooses a wavelength based upon the first fit
//					algorithm using quality
//
///////////////////////////////////////////////////////////////////
int ResourceManager::quality_first_fit_with_ordering(CreateConnectionProbeEvent* ccpe, unsigned int ci, bool* wave_available,size_t numberAvailableWaves)
{
	while(numberAvailableWaves > 0)
	{
		double xpm = 0.0;
		double fwm = 0.0;
		double ase = 0.0; 
		double Q_factor = 0.0;

		int wave = first_fit_with_ordering(ccpe,ci,wave_available);

		Q_factor = threadZero->getResourceManager()->estimate_Q(wave,ccpe->connectionPath,ccpe->connectionLength,&xpm,&fwm,&ase,ci);

		if(threadZero->getQualityParams().TH_Q <= Q_factor)
		{
			ccpe->wavelength = wave;

			print_connection_info(ccpe,Q_factor,ase,fwm,xpm,ci);

			delete[] wave_available;
			return wave;
		}
		else
		{
			--numberAvailableWaves;
			wave_available[wave] = false;
		}
	}

	delete[] wave_available;

	ccpe->wavelength = QUALITY_FAILURE;

	print_connection_info(ccpe,0.0,0.0,0.0,0.0,ci);

	return QUALITY_FAILURE;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	least_quality_fit
// Description:		Chooses a wavelength based upon the least quality
//					fit algorithm using quality
//
///////////////////////////////////////////////////////////////////
int ResourceManager::least_quality_fit(CreateConnectionProbeEvent* ccpe, unsigned int ci, bool* wave_available)
{
	int minQualityWave = -1;
	double minQualityQFactor = std::numeric_limits<double>::infinity();

	double minXPM = 0.0;
	double minFWM = 0.0;
	double minASE = 0.0;

	for(unsigned int w = 0; w < threadZero->getNumberOfWavelengths(); ++w)
	{
		if(wave_available[w] == true)
		{
			double qfactor = 0.0;

			double xpm = 0.0;
			double fwm = 0.0;
			double ase = 0.0;

			qfactor = threadZero->getResourceManager()->estimate_Q(w,ccpe->connectionPath,ccpe->connectionLength,&xpm,&fwm,&ase,ci);

			if(qfactor < minQualityQFactor && qfactor >= threadZero->getQualityParams().TH_Q)
			{
				minQualityWave = w;
				minQualityQFactor = qfactor;
				
				minXPM = xpm;
				minFWM = fwm;
				minASE = ase;
			}
		}
	}

	delete[] wave_available;

	if(minQualityWave == -1)
	{
		ccpe->wavelength = QUALITY_FAILURE;

		print_connection_info(ccpe,0.0,0.0,0.0,0.0,ci);
	}
	else
	{
		ccpe->wavelength = minQualityWave;

		print_connection_info(ccpe,minQualityQFactor,minASE,minFWM,minXPM,ci);
	}

	return ccpe->wavelength;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	most_quality_fit
// Description:		Chooses a wavelength based upon the least quality
//					fit algorithm using quality
//
///////////////////////////////////////////////////////////////////
int ResourceManager::most_quality_fit(CreateConnectionProbeEvent* ccpe, unsigned int ci, bool* wave_available)
{
	int maxQualityWave = -1;
	double maxQualityQFactor = 0.0;

	double maxXPM = 0.0;
	double maxFWM = 0.0;
	double maxASE = 0.0;

	for(unsigned int w = 0; w < threadZero->getNumberOfWavelengths(); ++w)
	{
		if(wave_available[w] == true)
		{
			double qfactor = 0.0;

			double xpm = 0.0;
			double fwm = 0.0;
			double ase = 0.0;

			qfactor = threadZero->getResourceManager()->estimate_Q(w,ccpe->connectionPath,ccpe->connectionLength,&xpm,&fwm,&ase,ci);

			if(qfactor > maxQualityQFactor && qfactor >= threadZero->getQualityParams().TH_Q)
			{
				maxQualityWave = w;
				maxQualityQFactor = qfactor;
				
				maxXPM = xpm;
				maxFWM = fwm;
				maxASE = ase;
			}
		}
	}

	delete[] wave_available;

	if(maxQualityWave == -1)
	{
		ccpe->wavelength = QUALITY_FAILURE;

		print_connection_info(ccpe,0.0,0.0,0.0,0.0,ci);
	}
	else
	{
		ccpe->wavelength = maxQualityWave;

		print_connection_info(ccpe,maxQualityQFactor,maxASE,maxFWM,maxXPM,ci);
	}

	return ccpe->wavelength;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	quality_random_fit
// Description:		Chooses a wavelength based upon the random
//					algorithm using quality
//
///////////////////////////////////////////////////////////////////
int ResourceManager::quality_random_fit(CreateConnectionProbeEvent* ccpe, unsigned int ci, bool* wave_available,size_t numberAvailableWaves)
{
	while(numberAvailableWaves > 0)
	{
		double xpm = 0.0;
		double fwm = 0.0;
		double ase = 0.0;
		double Q_factor = 0.0;

		int wave = random_fit(ccpe,ci,wave_available,numberAvailableWaves);

		Q_factor = threadZero->getResourceManager()->estimate_Q(wave,ccpe->connectionPath,ccpe->connectionLength,&xpm,&fwm,&ase,ci);

		if(threadZero->getQualityParams().TH_Q <=
		   threadZero->getResourceManager()->estimate_Q(wave,ccpe->connectionPath,ccpe->connectionLength,&xpm,&fwm,&ase,ci))
		{
			ccpe->wavelength = wave;

			print_connection_info(ccpe,Q_factor,ase,fwm,xpm,ci);

			delete[] wave_available;
			return wave;
		}
		else
		{
			--numberAvailableWaves;
			wave_available[wave] = false;
		}
	}

	delete[] wave_available;

	ccpe->wavelength = QUALITY_FAILURE;

	print_connection_info(ccpe,0.0,0.0,0.0,0.0,ci);

	return QUALITY_FAILURE;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	quality_most_used
// Description:		Chooses a wavelength based upon the random
//					algorithm using quality
//
///////////////////////////////////////////////////////////////////
int ResourceManager::quality_most_used(CreateConnectionProbeEvent* ccpe, unsigned int ci, bool* wave_available,size_t numberAvailableWaves)
{
	while(numberAvailableWaves > 0)
	{
		double xpm = 0.0;
		double fwm = 0.0;
		double ase = 0.0;
		double Q_factor = 0.0;

		int wave = most_used(ccpe,ci,wave_available);

		Q_factor = threadZero->getResourceManager()->estimate_Q(wave,ccpe->connectionPath,ccpe->connectionLength,&xpm,&fwm,&ase,ci);

		if(threadZero->getQualityParams().TH_Q <= Q_factor)
		{
			ccpe->wavelength = wave;

			print_connection_info(ccpe,Q_factor,ase,fwm,xpm,ci);

			delete[] wave_available;
			return wave;
		}
		else
		{
			--numberAvailableWaves;
			wave_available[wave] = false;
		}
	}

	delete[] wave_available;

	ccpe->wavelength = QUALITY_FAILURE;

	print_connection_info(ccpe,0.0,0.0,0.0,0.0,ci);

	return QUALITY_FAILURE;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	precompute_fwm_fs
// Description:		Precomputes the fwm_fs structures so they aren't
//					repeatedly computed as the simulation is running.
//
///////////////////////////////////////////////////////////////////
void ResourceManager::precompute_fwm_fs(std::vector<int> &fwm_nums)
{
	fwm_fs = new std::vector<double*>[threadZero->getNumberOfWavelengths()];
	inter_indecies = new std::vector<int*>[threadZero->getNumberOfWavelengths()];

	for(unsigned int w = 0; w < threadZero->getNumberOfWavelengths(); ++w)
	{
		fwm_fs->push_back(new double[threadZero->getNumberOfWavelengths()]);
		inter_indecies->push_back(new int[threadZero->getNumberOfWavelengths()]);

		fwm_nums.push_back(build_FWM_fs((*fwm_fs)[w],(*inter_indecies)[w],w));
	}
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	precompute_fwm_combinations
// Description:		Precomputes the combinations that cause FWM so
//					they aren't repeatedly computed as the simulation
//					is running.
//
///////////////////////////////////////////////////////////////////
void ResourceManager::precompute_fwm_combinations()
{
	std::vector<int> fwm_nums;

	precompute_fwm_fs(fwm_nums);

	fwm_combinations = new std::vector<int>[threadZero->getNumberOfWavelengths()];

	for(unsigned int w = 0; w < threadZero->getNumberOfWavelengths(); ++w)
	{
		wave_combines(sys_fs[w],(*fwm_fs)[w],fwm_nums[w],fwm_combinations[w]);
	}
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	initSPMatrix
// Description:		Initializes the SP matrix so that the paths
//					can be stored.
//
///////////////////////////////////////////////////////////////////
void ResourceManager::initSPMatrix()
{
	if(SP_paths != 0)
		return;

	SP_paths = new kShortestPathReturn*[threadZero->getNumberOfRouters() * threadZero->getNumberOfRouters()];

	for(unsigned int p = 0; p < threadZero->getNumberOfRouters() * threadZero->getNumberOfRouters(); ++p)
	{
		SP_paths[p] = 0;
	}
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	freeSPMatrix
// Description:		Free the SP matrix so that the memory can be
//					reused.
//
///////////////////////////////////////////////////////////////////
void ResourceManager::freeSPMatrix()
{
	if(SP_paths == 0)
		return;

	for(unsigned int p = 0; p < threadZero->getNumberOfRouters() * threadZero->getNumberOfRouters(); ++p)
	{
		if(SP_paths[p] != 0)
		{
			delete[] SP_paths[p]->pathcost;
			delete[] SP_paths[p]->pathinfo;
			delete[] SP_paths[p]->pathlen;

			delete SP_paths[p];
		}
	}

	delete[] SP_paths;

	SP_paths = 0;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	build_KSP_EdgeList
// Description:		Builds the KSP edge list structure and saves
//					it so that we don't have to keep recomputing it.
//
///////////////////////////////////////////////////////////////////
void ResourceManager::build_KSP_EdgeList()
{
	kSP_edgeList = new kShortestPathEdges[threadZero->getNumberOfEdges()];

	unsigned int num = 0;

	for(unsigned int a = 0; a < threadZero->getNumberOfRouters(); ++a)
	{
		for(unsigned int b = 0; b < threadZero->getNumberOfRouters(); ++b)
		{
			int edgeID = threadZero->getRouterAt(a)->isAdjacentTo(b);

			if(edgeID >= 0)
			{
				kSP_edgeList[num].src_node = a;
				kSP_edgeList[num].dest_node = b;
				kSP_edgeList[num].edge_cost = threadZero->getRouterAt(a)->
					getEdgeByDestination(b)->getNumberOfSpans();

				++num;
			}
		}
	}

	if(num > threadZero->getNumberOfEdges())
	{
		threadZero->recordEvent("ERROR: More edges than we expected.",true,0);
		exit(ERROR_TOO_MANY_EDGES);
	}
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	print_connection_info
// Description:		Prints the information about the connection
//
///////////////////////////////////////////////////////////////////
void ResourceManager::print_connection_info(CreateConnectionProbeEvent* ccpe, double Q_factor, double ase, double fwm, double xpm, unsigned int ci)
{
	std::string line;

	line.append("SETUP CONNECTION: Routers[" + std::to_string(ccpe->connectionLength) + "]:");

	for(unsigned int r = 0; r < ccpe->connectionLength; ++r)
	{
		line.append(std::to_string(ccpe->connectionPath[r]->getSourceIndex()));

		if(r < ccpe->connectionLength)
			line.append(",");
	}

	line.append(std::to_string(ccpe->connectionPath[ccpe->connectionLength-1]->getDestinationIndex()));

	line.append(" Wavelength = " + std::to_string(ccpe->wavelength));
	line.append(" Session = " + std::to_string(ccpe->session));
	line.append(" Sequence = " + std::to_string(ccpe->sequence));

	if(ccpe->wavelength >= 0)
	{
		std::ostringstream buffer;

		buffer << "Q-factor = " << std::setprecision(4) << Q_factor << ", XPM-noise = " << std::setprecision(4) << xpm << 
			", FWM-noise = " << std::setprecision(4) << fwm << ", ASE-noise = " << std::setprecision(4) << ase << std::endl;

		line.append(buffer.str());
	}
	else
	{
		line.append("\n");
	}

	threadZero->recordEvent(line,false,ci);

	if(ccpe->wavelength >= 0)
	{
		threads[ci]->getGlobalStats().aseNoiseTotal += ase;
		threads[ci]->getGlobalStats().fwmNoiseTotal += fwm;
		threads[ci]->getGlobalStats().xpmNoiseTotal += xpm;

		size_t spans = 0;

		for(unsigned int r = 0; r < ccpe->connectionLength; ++r)
		{
			spans += ccpe->connectionPath[r]->getNumberOfSpans();
		}
	}
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	calc_min_spans
// Description:		Prints the information about the connection
//
///////////////////////////////////////////////////////////////////
void ResourceManager::calc_min_spans()
{
	printf("Calculating router distances...");

	span_distance = new unsigned int[threadZero->getNumberOfRouters() * threadZero->getNumberOfRouters()];

	unsigned int maxMinDistance = 0;

	for(unsigned int r1 = 0; r1 < threadZero->getNumberOfRouters(); ++r1)
	{
		for(unsigned int r2 = 0; r2 < threadZero->getNumberOfRouters(); ++r2)
		{
			if(r1 != r2)
			{
				span_distance[r1 * threadZero->getNumberOfRouters() + r2] = calculate_span_distance(r1,r2);

				if(span_distance[r1 * threadZero->getNumberOfRouters() + r2] > maxMinDistance)
				{
					maxMinDistance = span_distance[r1 * threadZero->getNumberOfRouters() + r2];
				}
			}
			else
			{
				span_distance[r1 * threadZero->getNumberOfRouters() + r2] = 0;
			}
		}
	}

	threadZero->setMinDuration(maxMinDistance);

	threadZero->setQFactorMin(maxMinDistance);

	printf("done.\n");
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	generateWaveOrdering
// Description:		Generates the list of wave ordering to be used
//					by the first fit with ordering algorithm
//
///////////////////////////////////////////////////////////////////
void ResourceManager::generateWaveOrdering()
{
	wave_ordering = new int[threadZero->getNumberOfWavelengths()];
	wave_ordering[0] = 0;

	if(threadZero->getNumberOfWavelengths() > 1)
		wave_ordering[1] = threadZero->getNumberOfWavelengths() - 1;

	std::list< std::pair<int,int>* > wavesList;

	for(unsigned int w = 1; w < threadZero->getNumberOfWavelengths() - 1; ++w)
	{
		std::pair<int, int> *p1 = new std::pair<int, int>;

		p1->first = w;
		p1->second = 0;

		wavesList.push_back(p1);
	}

	while(wavesList.size() > 0)
	{
		//Calculate the min of the upper and lower bounds
		for(std::list< std::pair<int,int>* >::iterator iter1 = wavesList.begin(); iter1 != wavesList.end(); ++iter1)
		{
			(*iter1)->second = std::min(getLowerBound((*iter1)->first,static_cast<int>(threadZero->getNumberOfWavelengths()-wavesList.size())),
				getUpperBound((*iter1)->first,static_cast<int>(threadZero->getNumberOfWavelengths()-wavesList.size())));
		}

		std::list< std::pair<int,int>* >::iterator maxPair = wavesList.begin();

		//Find the wave with the maximum distance
		for(std::list< std::pair<int,int>* >::iterator iter2 = wavesList.begin(); iter2 != wavesList.end(); ++iter2)
		{
			if((*iter2)->second > (*maxPair)->second)
			{
				maxPair = iter2;
			}
		}

		wave_ordering[threadZero->getNumberOfWavelengths() - wavesList.size()] = (*maxPair)->first;

		delete *maxPair;

		wavesList.erase(maxPair);
	}
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	getLowerBound
// Description:		Calculates the wavelength on the left side of
//					w that has already been ordered
//
///////////////////////////////////////////////////////////////////
size_t ResourceManager::getLowerBound(int w, int n)
{
	size_t retval = threadZero->getNumberOfWavelengths();

	for(int a = 0; a < n; ++a)
	{
		if(wave_ordering[a] < w)
		{
			if(abs(wave_ordering[a] - w) < retval)
				retval = abs(wave_ordering[a] - w);
		}
	}

	return retval;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	getUpperBound
// Description:		Calculates the wavelength on the right side of
//					w that has already been ordered
//
///////////////////////////////////////////////////////////////////
size_t ResourceManager::getUpperBound(int w, int n)
{
	size_t retval = threadZero->getNumberOfWavelengths();

	for(int a = 0; a < n; ++a)
	{
		if(wave_ordering[a] > w)
		{
			if(abs(wave_ordering[a] - w) < retval)
				retval = abs(wave_ordering[a] - w);
		}
	}

	return retval;
}
