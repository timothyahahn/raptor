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

  for (size_t i = 0; i < threadZero->getNumberOfWavelengths(); i++)
    for (size_t j = 0; j < threadZero->getNumberOfWavelengths(); j++)
      sys_link_xpm_database[i * threadZero->getNumberOfWavelengths() + j] = 0.0;

#ifdef NO_OCTAVE
  threadZero->recordEvent(
      "Octave support not enabled, no xpm calculations possible", true, 0);
  return;
#endif  // NO_OCTAVE

  build_xpm_database(
      sys_fs, sys_fs_num, threadZero->getQualityParams().channel_power,
      threadZero->getQualityParams().D, threadZero->getQualityParams().alphaDB,
      threadZero->getQualityParams().gamma, res_disp);

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
	/*string_vector argvv(2);
	argvv(0) = "embedded";
	argvv(1) = "-q";
	octave_main(2, argvv.c_str_vec(), true);
	std::string warn_for;
	std::cerr << "about to call source_file" << std::endl;
	source_file("octave/hello_world.m",std::string(),true,true,warn_for);
	std::cerr << "returned from call source_file" << std::endl;
	std::cerr << "warn_for = " << warn_for << std::endl;
	octave_value_list inputs;
	const octave_value_list result = feval("hello_world", inputs, 0);
	result(0).print_raw(std::cout, true);
	clean_up_and_exit(0);*/
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
		octave_value_list inputs;
		const octave_value_list result = octave::feval("hello_world", inputs, 0);
		if (result.length() > 0)
		{
			std::cout << "hello_world returned" << result(0).int_value();
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
                                       double res_disp) {
  /*try {
          mwArray in1(1,fs_num,mxDOUBLE_CLASS);
          mwArray in2(1,1,mxDOUBLE_CLASS);
          mwArray in3(1,1,mxDOUBLE_CLASS);
          mwArray in4(1,1,mxDOUBLE_CLASS);
          mwArray in5(1,1,mxDOUBLE_CLASS);
          mwArray in6(1,1,mxDOUBLE_CLASS);
          mwArray in7(1,1,mxDOUBLE_CLASS);

          for(int id = 0;id<fs_num;id++)
                  in1(1,id+1) = fs[id];

          in2(1,1) = channel_power;
          in3(1,1) = D;
          in4(1,1) = alphaDB;
          in5(1,1) = gamma;
          in6(1,1) = res_disp;
          in7(1,1) = threadZero->getQualityParams().nonlinear_halfwin;

          int numberOfReturnVals = 1;
          mwArray out1(1,1,mxDOUBLE_CLASS);

          check_last_inputs(numberOfReturnVals,out1,in1,in2,in3,in4,in5,in6,in7);

          if(static_cast<int>(out1(1,1)) == 1)
          {
                  //The previous inputs are identical. Just reload the matrix.
                  printf("Reloading XPM matrix....");
          }
          else
          {
                  //The previous inputs were different. Rebuild the matrix.
                  printf("Building XPM matrix...");
                  build_libxpm_database(in1,in2,in3,in4,in5,in6,in7);
          }

  }
catch (const mwException& e) {
          std::cerr << e.what() << std::endl;
}
  catch (...) {
          std::cerr << "Unexpected error thrown" << std::endl;
}*/

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
  /*mwArray
  out(threadZero->getNumberOfWavelengths(),threadZero->getNumberOfWavelengths(),mxDOUBLE_CLASS);
  load_libxpm_database(1, out);

  for(int i=0;i<fs_num;i++)
          for(int j=0;j<fs_num;j++)
                  store[i * threadZero->getNumberOfWavelengths() + j] =
  out(i+1,j+1);

  printf("done.\n");*/

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
