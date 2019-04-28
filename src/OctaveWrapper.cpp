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

#ifndef NO_OCTAVE
#include <octave/oct.h>
#include <octave/octave.h>
#include <octave/parse.h>
#include <octave/interpreter.h>
#endif  // NO_OCTAVE

#include "Thread.h"

#include <algorithm>
#include <string>

extern Thread* threadZero;

double OctaveWrapper::res_disp = 0.0;

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
  return;
#endif  // NO_OCTAVE

  int identical = check_last_inputs(
	  sys_fs, sys_fs_num, threadZero->getQualityParams().channel_power,
	  threadZero->getQualityParams().D, threadZero->getQualityParams().alphaDB,
	  threadZero->getQualityParams().gamma, res_disp, threadZero->getQualityParams().halfwavelength);

  if (identical == 1)
  {
	  threadZero->recordEvent("XPM Matrix is identical, will not recompute.",true,0);
	  return;
  }
  else
  {
	  threadZero->recordEvent("XPM Matrix is diifferent, will recompute.", true, 0);
  }

  build_xpm_database(
      sys_fs, sys_fs_num, threadZero->getQualityParams().channel_power,
      threadZero->getQualityParams().D, threadZero->getQualityParams().alphaDB,
      threadZero->getQualityParams().gamma, res_disp, threadZero->getQualityParams().halfwavelength);

  load_xpm_database(sys_link_xpm_database, sys_fs_num);

  return;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	helloWorld
// Description:		Simple function to test octave
//
//
///////////////////////////////////////////////////////////////////
void OctaveWrapper::helloWorld()
{
#ifndef NO_OCTAVE
	try
	{
		octave::interpreter interp;
		interp.initialize();
		if (!interp.initialized())
		{
			std::cerr << "ERROR: Interpreter initialization failed" << std::endl;
			return;
		}
		int status = interp.execute();
		if (status != 0)
		{
			std::cerr << "ERROR: Creating embedded interpreter failed" << std::endl;
			return;
		}
		interp.get_load_path().append("octave",true);
		octave_value_list inputs;
		const octave_value_list result = octave::feval("hello_world", inputs, 0);
		if (result.length() > 0)
		{
			std::cout << "hello_world returned " << result(0).int_value() << std::endl;
		}
		else
		{
			std::cerr << "ERROR: result.length() is 0" << std::endl;
		}
	}
	catch (const octave::exit_exception & ex)
	{
		std::cerr << "Octave interpreter exited with status = "
			<< ex.exit_status() << std::endl;
	}
	catch (const octave::execution_exception&)
	{
		std::cerr << "error encountered in Octave evaluator!" << std::endl;
	}
#endif  // NO_OCTAVE
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
	try
	{
		octave::interpreter interp;
		interp.verbose(true);
		interp.initialize();
		if (!interp.initialized())
		{
			std::cerr << "ERROR: Interpreter initialization failed" << std::endl;
			return -1;
		}
		int status = interp.execute();
		if (status != 0)
		{
			std::cerr << "ERROR: Creating embedded interpreter failed" << std::endl;
			return -1;
		}
		interp.get_load_path().append("octave", true);
		octave_value_list inputs;
		NDArray fs_array(fs_num);
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
		const octave_value_list result = octave::feval("check_last_inputs", inputs, 7);
		if (result.length() > 0)
		{
			std::cout << "check_last_inputs returned " << result(0).int_value() << std::endl;
			return result(0).int_value();
		}
		else
		{
			std::cerr << "ERROR: result.length() is 0" << std::endl;
		}
	}
	catch (const octave::exit_exception & ex)
	{
		std::cerr << "Octave interpreter exited with status = "
			<< ex.exit_status() << std::endl;
	}
	catch (const octave::execution_exception&)
	{
		std::cerr << "error encountered in Octave evaluator!" << std::endl;
	}
#endif  // NO_OCTAVE
	return -1;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	build_xpm_database
// Description:		Builds a nonlinear database used to calculate
//					the XPM through a series of Octave
//calls.
//
///////////////////////////////////////////////////////////////////
void OctaveWrapper::build_xpm_database(double* fs, int fs_num,
                                       double channel_power, double D,
                                       double alphaDB, double gamma,
                                       double res_disp, double half_win) {
#ifndef NO_OCTAVE
	try
	{
		octave::interpreter interp;
		interp.verbose(true);
		interp.initialize();
		if (!interp.initialized())
		{
			std::cerr << "ERROR: Interpreter initialization failed" << std::endl;
			return;
		}
		int status = interp.execute();
		if (status != 0)
		{
			std::cerr << "ERROR: Creating embedded interpreter failed" << std::endl;
			return;
		}
		interp.get_load_path().append("octave", true);
		octave_value_list inputs;
		NDArray fs_array(fs_num);
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
		octave_value_list result = octave::feval("build_libxpm_database", inputs, 7);
	}
	catch (const octave::exit_exception & ex)
	{
		std::cerr << "Octave interpreter exited with status = "
			<< ex.exit_status() << std::endl;
	}
	catch (const octave::execution_exception&)
	{
		std::cerr << "error encountered in Octave evaluator!" << std::endl;
	}
#endif  // NO_OCTAVE
  return;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	load_xpm_database
// Description:		Loads the XPM database from the out matrix and
//					stores it into the store database.
//
///////////////////////////////////////////////////////////////////
void OctaveWrapper::load_xpm_database(double* store, int fs_num) {
#ifndef NO_OCTAVE
	try
	{
		octave::interpreter interp;
		interp.verbose(true);
		interp.initialize();
		if (!interp.initialized())
		{
			std::cerr << "ERROR: Interpreter initialization failed" << std::endl;
			return;
		}
		int status = interp.execute();
		if (status != 0)
		{
			std::cerr << "ERROR: Creating embedded interpreter failed" << std::endl;
			return;
		}
		interp.get_load_path().append("octave", true);
		octave_value_list inputs;
		octave_value_list result = octave::feval("load_xpm_database", inputs, 0);
	}
	catch (const octave::exit_exception & ex)
	{
		std::cerr << "Octave interpreter exited with status = "
			<< ex.exit_status() << std::endl;
	}
	catch (const octave::execution_exception&)
	{
		std::cerr << "error encountered in Octave evaluator!" << std::endl;
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
