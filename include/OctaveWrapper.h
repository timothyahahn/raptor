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

#ifndef NO_OCTAVE
#include <octave/oct.h>
#include <octave/octave.h>
#include <octave/parse.h>
#include <octave/interpreter.h>
#endif  // NO_OCTAVE

class OctaveWrapper {
 public:
  OctaveWrapper();

  void build_nonlinear_datastructure(double *sys_fs,
                                            double *sys_link_xpm_database);
  void helloWorld();

 private:
  void build_xpm_database(double *fs, int fs_num, double channel_power,
                                 double D, double alphaDB, double gamma,
                                 double res_disp, double half_win);
  int check_last_inputs(double* fs, int fs_num, double channel_power,
	                           double D, double alphaDB, double gamma,
	                           double res_disp, double half_win);
  void load_xpm_database(double *store, int fs_num);
  int gen_frequency_comb(double *frequencies, double fc, double step,
                                int left, int right, int wo_fc);

#ifndef NO_OCTAVE
  octave::interpreter interp;
#endif //NO_OCTAVE

  double res_disp;
};

#endif
