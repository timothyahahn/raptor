// ____________________________________________________________________________
//
//  General Information:
//
//  File Name:      OctaveWrapper.cpp
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

#include "OctaveWrapper.h"

#include "Thread.h"

#include <algorithm>
#include <string>

extern Thread* threadZero;

OctaveWrapper::OctaveWrapper()
	: res_disp(0.0)
{
#ifndef NO_OCTAVE
	threadZero->recordEvent("Initializing Octave Interpreter", true, 0);
	threadZero->flushLog(true);
	try
	{
		interp.initialize();
		if (!interp.initialized())
		{
			threadZero->recordEvent("ERROR: Interpreter initialization failed", true, 0);
			threadZero->flushLog(true);
			return;
		}
		int status = interp.execute();
		if (status != 0)
		{
			threadZero->recordEvent("ERROR: Creating embedded interpreter failed", true, 0);
			threadZero->flushLog(true);
			return;
		}
		interp.get_load_path().append("octave", true);
	}
	catch (const octave::exit_exception & ex)
	{
		threadZero->recordEvent("ERROR: Octave interpreter exited with status = " + ex.exit_status(), true, 0);
		threadZero->flushLog(true);
	}
	catch (const octave::execution_exception&)
	{
		threadZero->recordEvent("ERROR: encountered in Octave evaluator!", true, 0);
		threadZero->flushLog(true);
	}
	catch (std::runtime_error & re)
	{
		threadZero->recordEvent("ERROR: " + std::string(re.what()), true, 0);
		threadZero->flushLog(true);
	}
#endif  // NO_OCTAVE
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	build_nonlinear_datastructure
// Description:		Builds a nonlinear database used to calculate
//					the XPM through a series of Octave
//calls.
//
///////////////////////////////////////////////////////////////////
void OctaveWrapper::build_nonlinear_datastructure(
    double* sys_fs, double* sys_link_xpm_database) {
  int sys_fs_num =
      gen_frequency_comb(sys_fs, threadZero->getQualityParams().fc,
                         threadZero->getQualityParams().f_step,
                         threadZero->getQualityParams().halfwavelength,
                         threadZero->getQualityParams().halfwavelength, 1);

  std::fill(sys_link_xpm_database, sys_link_xpm_database + 
	  threadZero->getNumberOfWavelengths() * threadZero->getNumberOfWavelengths(), 0.0);

#ifdef NO_OCTAVE
  threadZero->recordEvent(
      "Octave support not enabled, no xpm calculations possible", true, 0);
  threadZero->flushLog(true);
  return;
#endif  // NO_OCTAVE

  int identical = check_last_inputs(
	  sys_fs, sys_fs_num, threadZero->getQualityParams().channel_power,
	  threadZero->getQualityParams().D, threadZero->getQualityParams().alphaDB,
	  threadZero->getQualityParams().gamma, res_disp, threadZero->getQualityParams().halfwavelength);

  if (identical == 1)
  {
	  threadZero->recordEvent("XPM Matrix is identical, will not recompute.",true,0);
	  threadZero->flushLog(true);
	  return;
  }
  else
  {
	  threadZero->recordEvent("XPM Matrix is different, will recompute.", true, 0);
	  threadZero->flushLog(true);
  }

  build_libxpm_database(
      sys_fs, sys_fs_num, threadZero->getQualityParams().channel_power,
      threadZero->getQualityParams().D, threadZero->getQualityParams().alphaDB,
      threadZero->getQualityParams().gamma, res_disp, threadZero->getQualityParams().halfwavelength);

  load_libxpm_database(sys_link_xpm_database, sys_fs_num);

  return;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	check_last_inputs
// Description:		Checks the nonlinear datastructure to determine
//                  if the inputs have changed and a recalculation
//                  is necessary.
//
///////////////////////////////////////////////////////////////////
int OctaveWrapper::check_last_inputs(double* fs, int fs_num,
	double channel_power, double D,
	double alphaDB, double gamma,
	double res_disp, double half_win) {
#ifndef NO_OCTAVE
	threadZero->recordEvent("Calling check_last_inputs.m", true, 0);
	threadZero->flushLog(true);
	try
	{
		if (!interp.initialized())
		{
			threadZero->recordEvent("ERROR: Interpreter not initialized", true, 0);
			threadZero->flushLog(true);
			return ERROR_OCTAVE;
		}
		octave_value_list inputs(7);
		dim_vector dv(fs_num,1);
		Matrix fs_array(dv);
		for (size_t f = 0; f < fs_num; ++f)
		{
			fs_array(f) = fs[f];
		}
		inputs(0) = fs_array;
		inputs(1) = channel_power;
		inputs(2) = D;
		inputs(3) = alphaDB;
		inputs(4) = gamma;
		inputs(5) = res_disp;
		inputs(6) = half_win;
		const octave_value_list result = octave::feval("check_last_inputs", inputs, 1);
		if (result.length() > 0)
		{
			return result(0).int_value();
		}
		else
		{
			threadZero->recordEvent("ERROR: result.length() is 0", true, 0);
			threadZero->flushLog(true);
		}
	}
	catch (const octave::exit_exception & ex)
	{
		threadZero->recordEvent("ERROR: Octave interpreter exited with status = " + ex.exit_status(), true, 0);
		threadZero->flushLog(true);
	}
	catch (const octave::execution_exception&)
	{
		threadZero->recordEvent("ERROR: encountered in Octave evaluator!", true, 0);
		threadZero->flushLog(true);
	}
	catch (std::runtime_error & re)
	{
		threadZero->recordEvent("ERROR: " + std::string(re.what()), true, 0);
		threadZero->flushLog(true);
	}
#endif  // NO_OCTAVE
	return ERROR_OCTAVE;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	build_libxpm_database
// Description:		Builds a nonlinear database used to calculate
//					the XPM through a series of Octave
//calls.
//
///////////////////////////////////////////////////////////////////
void OctaveWrapper::build_libxpm_database(double* fs, int fs_num,
                                       double channel_power, double D,
                                       double alphaDB, double gamma,
                                       double res_disp, double half_win) {
#ifndef NO_OCTAVE
	threadZero->recordEvent("Calling build_libxpm_database.m", true, 0);
	threadZero->flushLog(true);
	try
	{
		if (!interp.initialized())
		{
			threadZero->recordEvent("ERROR: Interpreter not initialized", true, 0);
			threadZero->flushLog(true);
			return;
		}
		octave_value_list inputs(7);
		dim_vector dv(fs_num, 1);
		Matrix fs_array(dv);
		for (size_t f = 0; f < fs_num; ++f)
		{
			fs_array(f) = fs[f];
		}
		inputs(0) = fs_array;
		inputs(1) = channel_power;
		inputs(2) = D;
		inputs(3) = alphaDB;
		inputs(4) = gamma;
		inputs(5) = res_disp;
		inputs(6) = half_win;
		octave_value_list result = octave::feval("build_libxpm_database", inputs, 0);
	}
	catch (const octave::exit_exception &ex)
	{
		threadZero->recordEvent("ERROR: Octave interpreter exited with status = " + ex.exit_status(), true, 0);
		threadZero->flushLog(true);
	}
	catch (const octave::execution_exception&)
	{
		threadZero->recordEvent("ERROR: encountered in Octave evaluator!", true, 0);
		threadZero->flushLog(true);
	}
	catch (std::runtime_error &re)
	{
		threadZero->recordEvent("ERROR: " + std::string(re.what()), true, 0);
		threadZero->flushLog(true);
	}
#endif  // NO_OCTAVE
  return;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	load_libxpm_database
// Description:		Loads the XPM database from the out matrix and
//					stores it into the store database.
//
///////////////////////////////////////////////////////////////////
void OctaveWrapper::load_libxpm_database(double* store, int fs_num) {
#ifndef NO_OCTAVE
	threadZero->recordEvent("Calling load_libxpm_database.m", true, 0);
	threadZero->flushLog(true);
	try
	{
		if (!interp.initialized())
		{
			threadZero->recordEvent("ERROR: Interpreter not initialized", true, 0);
			threadZero->flushLog(true);
			return;
		}
		octave_value_list result = octave::feval("load_libxpm_database");
		if (result.length() > 0)
		{
			Matrix m = result(0).array_value();

			for (int a = 0; a < fs_num; ++a)
				for (int b = 0; b < fs_num; ++b)
					store[a * fs_num + b] = m(a, b);
		}
		else
		{
			threadZero->recordEvent("ERROR: result.length() is unexpectedly 0", true, 0);
			threadZero->flushLog(true);
		}
	}
	catch (const octave::exit_exception &ex)
	{
		threadZero->recordEvent("ERROR: Octave interpreter exited with status = " + ex.exit_status(), true, 0);
		threadZero->flushLog(true);
	}
	catch (const octave::execution_exception&)
	{
		threadZero->recordEvent("ERROR: encountered in Octave evaluator!", true, 0);
		threadZero->flushLog(true);
	}
	catch (std::runtime_error &re)
	{
		threadZero->recordEvent("ERROR: " + std::string(re.what()), true, 0);
		threadZero->flushLog(true);
	}
#endif  // NO_OCTAVE
	return;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	gen_frequency_comb
// Description:		Generates the frequency combinations and returns
//					the number.
//
///////////////////////////////////////////////////////////////////
int OctaveWrapper::gen_frequency_comb(double* frequencies, double fc,
                                      double step, int left, int right,
                                      int wo_fc) {
  int i, num = 0;

  for (i = -left; i < 0; i++) frequencies[num++] = fc + i * step;

  if (wo_fc) frequencies[num++] = fc;

  for (i = 1; i <= right; i++) frequencies[num++] = fc + i * step;

  return num;
}
