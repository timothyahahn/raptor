// ____________________________________________________________________________
//
//  General Information:
//
//  File Name:      ResourceManager.h
//  Author:         Timothy Hahn, PhD
//  Project:        raptor
//
//  Description:    The file contains the declaration of the ResourceManager class.
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

#ifndef RESOURCE_MANAGER_H
#define RESOURCE_MANAGER_H

#include <queue>
#include <utility>
#include <vector>

#include "Edge.h"
#include "Event.h"
#include "Router.h"

#include "QYInclude.h"

using std::less;
using std::queue;
using std::pair;
using std::priority_queue;
using std::vector;

struct DP_item
{
	Edge** path;
	unsigned int pathLength;
	unsigned int pathSpans;
	bool* waveAvailability;
};

class ResourceManager
{
	public:
		ResourceManager();
		~ResourceManager();

		kShortestPathReturn* calculate_SP_path(unsigned int src, unsigned int dest, unsigned int k, unsigned int ci);
		kShortestPathReturn* calculate_LORA_path(unsigned int src, unsigned int dest, unsigned int k, unsigned int ci);
		kShortestPathReturn* calculate_PAR_path(unsigned int src_index, unsigned int dest_index, unsigned int k, unsigned int ci);
		kShortestPathReturn* calculate_IA_path(unsigned int src_index, unsigned int dest_index, unsigned int ci);
		kShortestPathReturn* calculate_QM_path(unsigned int src_index, unsigned int dest_index, unsigned int k, unsigned int ci);
		kShortestPathReturn* calculate_AQoS_path(unsigned int src_index, unsigned int dest_index, unsigned int k, unsigned int ci);
		kShortestPathReturn* calculate_DP_path(unsigned int src_index, unsigned int dest_index, unsigned int k, unsigned int ci);
		kShortestPathReturn* calculate_ACO_path(unsigned int src_index, unsigned int dest_index, unsigned int k, unsigned int ci);
		kShortestPathReturn* calculate_MM_ACO_path(unsigned int src_index, unsigned int dest_index, unsigned int k, unsigned int ci);

		int choose_wavelength(CreateConnectionProbeEvent* ccpe, unsigned int ci);

		double estimate_Q(int lambda, Edge **Path, unsigned int pathLen, double *xpm, double *fwm, double *ase, unsigned int ci);

		void initSPMatrix();
		void freeSPMatrix();

		double path_fwm_term(int spans,double fi,double fj, double fk,double fc,int dgen);
		double path_xpm_term(int spans, int lambda, int wave);

		void print_connection_info(CreateConnectionProbeEvent* ccpe, double Q_factor, double ase, double fwm, double xpm, unsigned int ci);

		double* sys_fs;
		vector<int>* fwm_combinations;

		unsigned int* span_distance;

	private:
		double path_ase_noise(int lambda, Edge **Path, unsigned int pathLen, unsigned int ci);

		double path_fwm_noise(int lambda, Edge **Path, unsigned int pathLen, unsigned int ci);

		double path_xpm_noise(int lambda, Edge **Path, unsigned int pathLen, unsigned int ci);
		
		int build_FWM_fs(double *inter_fs,int *inter_indecies, int lambda);
		int wave_combines(double fc, double *fs,int fs_num, vector<int> &fs_coms);
		bool can_find(int fi,int fj,int fk,vector<int> &fs_coms,int com_num);
		int degeneracy(int fi,int fj,int fk);

		double* sys_link_xpm_database;
		int sys_fs_num;

		int first_fit(CreateConnectionProbeEvent* ccpe, unsigned int ci, bool* wave_available);
		int first_fit_with_ordering(CreateConnectionProbeEvent* ccpe, unsigned int ci, bool* wave_available);

		int random_fit(CreateConnectionProbeEvent* ccpe, unsigned int ci, bool* wave_available,unsigned int numberAvailableWaves);

		int most_used(CreateConnectionProbeEvent* ccpe, unsigned int ci, bool* wave_available);

		int quality_first_fit(CreateConnectionProbeEvent* ccpe, unsigned int ci, bool* wave_available,unsigned int numberAvailableWaves);
		int quality_first_fit_with_ordering(CreateConnectionProbeEvent* ccpe, unsigned int ci, bool* wave_available,unsigned int numberAvailableWaves);

		int quality_random_fit(CreateConnectionProbeEvent* ccpe, unsigned int ci, bool* wave_available,unsigned int numberAvailableWaves);

		int quality_most_used(CreateConnectionProbeEvent* ccpe, unsigned int ci, bool* wave_available,unsigned int numberAvailableWaves);

		int least_quality_fit(CreateConnectionProbeEvent* ccpe, unsigned int ci, bool* wave_available);
		int most_quality_fit(CreateConnectionProbeEvent* ccpe, unsigned int ci, bool* wave_available);

		void precompute_fwm_fs(vector<int> &fwm_nums);
		void precompute_fwm_combinations();

		vector <double*>* fwm_fs;
		vector <int*>* inter_indecies;

		kShortestPathReturn** SP_paths;

		void build_KSP_EdgeList();

		kShortestPathEdges* kSP_edgeList;

		void calc_min_spans();
		unsigned int calculate_span_distance(unsigned int src, unsigned int dest);

		int* wave_ordering;

		void generateWaveOrdering();

		int getLowerBound(int w, int n);
		int getUpperBound(int w, int n);
};

struct Ant
{
	Router* location;
	Edge** path;
	int pathlen;
	double Q;
};

#endif
