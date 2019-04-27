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

struct DP_item
{
	Edge** path;
	unsigned int pathLength;
	size_t pathSpans;
	bool* waveAvailability;
};

class ResourceManager
{
	public:
		ResourceManager();
		~ResourceManager();

		kShortestPathReturn* calculate_SP_path(size_t src_index, size_t dest_index, size_t k, size_t ci);
		kShortestPathReturn* calculate_LORA_path(size_t src_index, size_t dest_index, size_t k, size_t ci);
		kShortestPathReturn* calculate_PAR_path(size_t src_index, size_t dest_index, size_t k, size_t ci);
		kShortestPathReturn* calculate_IA_path(size_t src_index, size_t dest_index, size_t ci);
		kShortestPathReturn* calculate_QM_path(size_t src_index, size_t dest_index, size_t k, size_t ci);
		kShortestPathReturn* calculate_AQoS_path(size_t src_index, size_t dest_index, size_t k, size_t ci);
		kShortestPathReturn* calculate_DP_path(size_t src_index, size_t dest_index, size_t k, size_t ci);
		kShortestPathReturn* calculate_ACO_path(size_t src_index, size_t dest_index, size_t k, size_t ci);
		kShortestPathReturn* calculate_MM_ACO_path(size_t src_index, size_t dest_index, size_t k, size_t ci);

		long long int choose_wavelength(CreateConnectionProbeEvent* ccpe, size_t ci);

		double estimate_Q(long long int lambda, Edge **Path, size_t pathLen, double *xpm, double *fwm, double *ase, size_t ci);

		void initSPMatrix();
		void freeSPMatrix();

		double path_fwm_term(size_t spans,double fi,double fj, double fk,double fc,int dgen);
		double path_xpm_term(size_t spans, size_t lambda, size_t wave);

		void print_connection_info(CreateConnectionProbeEvent* ccpe, double Q_factor, double ase, double fwm, double xpm, size_t ci);

		double* sys_fs;
		std::vector<int>* fwm_combinations;

		size_t* span_distance;

	private:
		double path_ase_noise(long long int lambda, Edge **Path, size_t pathLen, size_t ci);

		double path_fwm_noise(long long int lambda, Edge **Path, size_t pathLen, size_t ci);

		double path_xpm_noise(long long int lambda, Edge **Path, size_t pathLen, size_t ci);
		
		int build_FWM_fs(double *inter_fs,int *inter_indecies, int lambda);
		int wave_combines(double fc, double *fs,int fs_num, std::vector<int> &fs_coms);
		bool can_find(int fi,int fj,int fk,std::vector<int> &fs_coms,int com_num);
		int degeneracy(int fi,int fj,int fk);

		double* sys_link_xpm_database;
		int sys_fs_num;

		long long int first_fit(CreateConnectionProbeEvent* ccpe, size_t ci, bool* wave_available);
		long long int first_fit_with_ordering(CreateConnectionProbeEvent* ccpe, size_t ci, bool* wave_available);

		long long int random_fit(CreateConnectionProbeEvent* ccpe, size_t ci, bool* wave_available, size_t numberAvailableWaves);

		long long int most_used(CreateConnectionProbeEvent* ccpe, size_t ci, bool* wave_available);

		long long int quality_first_fit(CreateConnectionProbeEvent* ccpe, size_t ci, bool* wave_available,size_t numberAvailableWaves);
		long long int quality_first_fit_with_ordering(CreateConnectionProbeEvent* ccpe, size_t ci, bool* wave_available,size_t numberAvailableWaves);

		long long int quality_random_fit(CreateConnectionProbeEvent* ccpe, size_t ci, bool* wave_available,size_t numberAvailableWaves);

		long long int quality_most_used(CreateConnectionProbeEvent* ccpe, size_t ci, bool* wave_available,size_t numberAvailableWaves);

		long long int least_quality_fit(CreateConnectionProbeEvent* ccpe, size_t ci, bool* wave_available);
		long long int most_quality_fit(CreateConnectionProbeEvent* ccpe, size_t ci, bool* wave_available);

		void precompute_fwm_fs(std::vector<int> &fwm_nums);
		void precompute_fwm_combinations();

		std::vector <double*>* fwm_fs;
		std::vector <int*>* inter_indecies;

		kShortestPathReturn** SP_paths;

		void build_KSP_EdgeList();

		kShortestPathEdges* kSP_edgeList;

		void calc_min_spans();
		unsigned int calculate_span_distance(size_t src_index, size_t dest_index);

		long long int* wave_ordering;

		void generateWaveOrdering();

		long long int getLowerBound(int w, int n);
		long long int getUpperBound(int w, int n);
};

struct Ant
{
	Router* location;
	Edge** path;
	int pathlen;
	double Q;
};

#endif
