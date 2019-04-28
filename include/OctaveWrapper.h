// ____________________________________________________________________________
//
//  General Information:
//
//  File Name:      OctaveWrapper.h
//  Author:         Timothy Hahn, PhD
//  Project:        raptor
//
//  Description:    The file contains the calls to Octave, required to calculate
//  the
//					XPM modules.
//
//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//  Revision History:
//
//  04/14/2019  v2.0    Reworked version based upon cmake and octave
//
// ____________________________________________________________________________

#ifndef OCTAVE_WRAPPER_H
#define OCTAVE_WRAPPER_H

class OctaveWrapper {
 public:
  static void build_nonlinear_datastructure(double *sys_fs,
                                            double *sys_link_xpm_database);
  static void helloWorld();

 private:
  static void build_xpm_database(double *fs, int fs_num, double channel_power,
                                 double D, double alphaDB, double gamma,
                                 double res_disp, double half_win);
  static int check_last_inputs(double* fs, int fs_num, double channel_power,
	                           double D, double alphaDB, double gamma,
	                           double res_disp, double half_win);
  static void load_xpm_database(double *store, int fs_num);
  static int gen_frequency_comb(double *frequencies, double fc, double step,
                                int left, int right, int wo_fc);

  static double res_disp;
};

#endif
