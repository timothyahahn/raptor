// ____________________________________________________________________________
//
//  General Information:
//
//  File Name:      ResourceManager.h
//  Author:         Timothy Hahn, PhD
//  Project:        raptor
//
//  Description:    The file contains the declaration of the ResourceManager
//  class.
//					The purpose of the ResourceManager is to
//calculate the path 					from source to destination, calculate the Q-factor, and
//					track the availability of the network
//resources.
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

struct DP_item {
  Edge** path;
  size_t pathLength;
  size_t pathSpans;
  bool* waveAvailability;
};

class ResourceManager {
 public:
  ResourceManager();
  ~ResourceManager();

  kShortestPathReturn* calculate_SP_path(size_t src_index, size_t dest_index,
                                         size_t k, size_t ci);
  kShortestPathReturn* calculate_LORA_path(size_t src_index, size_t dest_index,
                                           size_t k, size_t ci);
  kShortestPathReturn* calculate_PAR_path(size_t src_index, size_t dest_index,
                                          size_t k, size_t ci);
  kShortestPathReturn* calculate_IA_path(size_t src_index, size_t dest_index,
                                         size_t ci);
  kShortestPathReturn* calculate_QM_path(size_t src_index, size_t dest_index,
                                         size_t k, size_t ci);
  kShortestPathReturn* calculate_AQoS_path(size_t src_index, size_t dest_index,
                                           size_t k, size_t ci);
  kShortestPathReturn* calculate_DP_path(size_t src_index, size_t dest_index,
                                         size_t k, size_t ci);
  kShortestPathReturn* calculate_ACO_path(size_t src_index, size_t dest_index,
                                          size_t k, size_t ci);
  kShortestPathReturn* calculate_MM_ACO_path(size_t src_index,
                                             size_t dest_index, size_t k,
                                             size_t ci);

  long long int choose_wavelength(CreateConnectionProbeEvent* ccpe, size_t ci);

  double estimate_Q(long long int lambda, Edge** Path, size_t pathLen,
                    double* xpm, double* fwm, double* ase, size_t ci) const;

  void initSPMatrix();
  void freeSPMatrix();

  void print_connection_info(CreateConnectionProbeEvent* ccpe, double Q_factor,
                             double ase, double fwm, double xpm, size_t ci) const;

  size_t* span_distance;

 private:
  double path_ase_noise(long long int lambda, Edge** Path, size_t pathLen,
                        size_t ci) const;

  double path_fwm_noise(long long int lambda, Edge** Path, size_t pathLen,
                        size_t ci) const;

  double path_xpm_noise(long long int lambda, Edge** Path, size_t pathLen,
                        size_t ci) const;

  double path_fwm_term(size_t spans, double fi, double fj, double fk, double fc,
	  long long int dgen) const;
  double path_xpm_term(size_t spans, size_t lambda, size_t wave) const;

  size_t calculate_span_distance(size_t src_index, size_t dest_index);

  long long int build_FWM_fs(double* inter_fs, long long int* inter_indecies, size_t lambda);
  long long int wave_combines(double fc, double* fs, long long int fs_num,
                    std::vector<long long int>& fs_coms);
  bool can_find(long long int fi, long long int fj, long long int fk, std::vector<long long int>& fs_coms, long long int com_num);
  long long int degeneracy(long long int fi, long long int fj, long long int fk);

  double* sys_link_xpm_database;
  long long int sys_fs_num;

  long long int first_fit(CreateConnectionProbeEvent* ccpe, size_t ci,
                          bool* wave_available);
  long long int first_fit_with_ordering(CreateConnectionProbeEvent* ccpe,
                                        size_t ci, bool* wave_available);

  long long int random_fit(CreateConnectionProbeEvent* ccpe, size_t ci,
                           bool* wave_available, size_t numberAvailableWaves);

  long long int most_used(CreateConnectionProbeEvent* ccpe, size_t ci,
                          bool* wave_available);

  long long int quality_first_fit(CreateConnectionProbeEvent* ccpe, size_t ci,
                                  bool* wave_available,
                                  size_t numberAvailableWaves);
  long long int quality_first_fit_with_ordering(
      CreateConnectionProbeEvent* ccpe, size_t ci, bool* wave_available,
      size_t numberAvailableWaves);

  long long int quality_random_fit(CreateConnectionProbeEvent* ccpe, size_t ci,
                                   bool* wave_available,
                                   size_t numberAvailableWaves);

  long long int quality_most_used(CreateConnectionProbeEvent* ccpe, size_t ci,
                                  bool* wave_available,
                                  size_t numberAvailableWaves);

  long long int least_quality_fit(CreateConnectionProbeEvent* ccpe, size_t ci,
                                  bool* wave_available);
  long long int most_quality_fit(CreateConnectionProbeEvent* ccpe, size_t ci,
                                 bool* wave_available);

  void precompute_fwm_fs(std::vector<long long int>& fwm_nums);
  void precompute_fwm_combinations();

  std::vector<double*>* fwm_fs;
  std::vector<long long int*>* inter_indecies;

  kShortestPathReturn** SP_paths;

  void build_KSP_EdgeList();

  kShortestPathEdges* kSP_edgeList;

  void calc_min_spans();

  long long int* wave_ordering;

  void generateWaveOrdering();

  double* sys_fs;
  std::vector<long long int>* fwm_combinations;

  long long int getLowerBound(int w, int n);
  long long int getUpperBound(int w, int n);
};

struct Ant {
  Router* location;
  Edge** path;
  int pathlen;
  double Q;
};

#endif
