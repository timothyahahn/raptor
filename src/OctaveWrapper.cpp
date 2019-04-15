// ____________________________________________________________________________
//
//  General Information:
//
//  File Name:      OctaveWrapper.cpp
//  Author:         Timothy Hahn, PhD
//  Project:        raptor
//
//  Description:    The file contains the calls to Octave, required to calculate the
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

extern Thread* threadZero;

double OctaveWrapper::res_disp = 0.0;

///////////////////////////////////////////////////////////////////
//
// Function Name:	build_nonlinear_datastructure
// Description:		Builds a nonlinear database used to calculate
//					the XPM through a series of Octave calls.
//
///////////////////////////////////////////////////////////////////
void OctaveWrapper::build_nonlinear_datastructure(double* sys_fs, double* sys_link_xpm_database)
{
	int sys_fs_num = gen_frequency_comb(sys_fs,threadZero->getQualityParams().fc,
		threadZero->getQualityParams().f_step,threadZero->getQualityParams().halfwavelength,
		threadZero->getQualityParams().halfwavelength,1);

	for(unsigned short int i = 0; i < threadZero->getNumberOfWavelengths(); i++)
		for(unsigned short int j = 0; j < threadZero->getNumberOfWavelengths(); j++)
			sys_link_xpm_database[i * threadZero->getNumberOfWavelengths() + j] = 0.0;

	build_xpm_database(sys_fs,sys_fs_num,threadZero->getQualityParams().channel_power,
		threadZero->getQualityParams().D,threadZero->getQualityParams().alphaDB,
		threadZero->getQualityParams().gamma,res_disp);

	load_xpm_database(sys_link_xpm_database,sys_fs_num);

	return;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	build_xpm_database
// Description:		Builds a nonlinear database used to calculate
//					the XPM through a series of Octave calls.
//
///////////////////////////////////////////////////////////////////
void OctaveWrapper::build_xpm_database(double *fs, int fs_num,double channel_power,double D,double alphaDB,double gamma,double res_disp)
{
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
void OctaveWrapper::load_xpm_database(double* store,int fs_num)
{
	/*mwArray out(threadZero->getNumberOfWavelengths(),threadZero->getNumberOfWavelengths(),mxDOUBLE_CLASS);
	load_libxpm_database(1, out);

	for(int i=0;i<fs_num;i++)
		for(int j=0;j<fs_num;j++)
			store[i * threadZero->getNumberOfWavelengths() + j] = out(i+1,j+1);

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
int OctaveWrapper::gen_frequency_comb(double *frequencies,double fc,double step, int left,int right, int wo_fc)
{
	int i, num = 0;

	for(i=-left;i<0;i++)
		frequencies[num++]= fc+i*step;

	if (wo_fc)
		frequencies[num++] = fc;

	for(i=1;i<= right;i++)
		frequencies[num++]= fc+i*step;

	return num;
}
