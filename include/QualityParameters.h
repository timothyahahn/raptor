// ____________________________________________________________________________
//
//  General Information:
//
//  File Name:      QualityParameters.h
//  Author:         Timothy Hahn, PhD
//  Project:        raptor
//
//  Description:    The file contains the declaration of the QualityParameters
//					struct, which are used to calculate the Q
//factor
//
//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//  Revision History:
//
//  05/20/2009	v1.0	Initial Version.
//  04/14/2019  v2.0    Reworked version based upon cmake and octave
//
// ____________________________________________________________________________

#ifndef QUALITY_PARAMS_H
#define QUALITY_PARAMS_H

enum DestinationDistribution {
  UNIFORM = 1,
  DISTANCE = 2,
  INVERSE_DISTANCE = 3
};

struct QualityParameters {
  double arrival_interval;       // the inter arrival time on each workstation
  double duration;               // the duration time of one connection
  int nonlinear_halfwin;         // the total nonlinear_win is 2*M+1
  int halfwavelength;            // the total wavelength is 2N+1
  double fc;                     // center of the wavelength comb
  double f_step;                 // step of the wavelength comb
  double channel_power;          // power per channel
  double L;                      // length of the NZDSF in each span
  double alphaDB;                // attenuation of NZDSF in each span
  double alpha;                  // attenuation of NZDSF in each span
  double D;                      // dispersion of NZDSF in each span
  double S;                      // disperison slop of NZDSF in each span
  double gamma;                  // nonlinear coefficent of NZDSF in each span
  double QFactor_factor;         // Factor of max min to set the threshold
  double TH_Q;                   // nonlinear coefficent of NZDSF in each span
  double EDFA_Noise_Figure;      // EDFA Noise Figure
  double EDFA_Gain;              // EDFA Gain
  double B_w;                    // Optical bandwidth of the signal channel
  double *ASE_perEDFA;           // ASE noise per EDFA
  double usage_update_interval;  // interval for updating usage for PABR and
                                 // LORA
  double beta;              // beta value for PABR and LORA
  int gui_update_interval;  // interval for updating the gui
  unsigned int max_probes;  // max number of probes per connection request
  double refractive_index;  // refractive index of the optical links
  bool q_factor_stats;      // should the program calculate the Q-factor stats
                        // (1=yes,0=no)
  bool detailed_log;  // should the program keep a detailed log (1=yes,0=no)
  DestinationDistribution dest_dist;  // distribution of the destination
  double DP_alpha;                    // Alpha value for Dynamic Programming
  unsigned int ACO_ants;              // number of ants in each ACO iteration
  double ACO_alpha;                   // the pheromone power index for ACO
  double ACO_beta;      // the heuristic information power index for ACO
  double ACO_rho;       // the pheromone evaporation rate for ACO
  double MM_ACO_gamma;  // the min-max pheromone ratio for MM ACO
  unsigned int
      MM_ACO_N_iter;  // the number of iterations for stagnation for MM ACO
  unsigned int
      MM_ACO_N_reset;  // the number of reinitialization times for MM ACO
};

#endif
