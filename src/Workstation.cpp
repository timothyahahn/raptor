// ____________________________________________________________________________
//
//  General Information:
//
//  File Name:      Workstation.cpp
//  Author:         Timothy Hahn, PhD
//  Project:        raptor
//
//  Description:    The file contains the implementation of the Workstation class.
//					The purpose of the Workstations are to generate traffic for
//					the optical network.
//
//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//  Revision History:
//
//  05/20/2009	v1.0	Initial Version.
//  04/14/2019  v2.0    Reworked version based upon cmake and octave
//
// ____________________________________________________________________________

#include "Workstation.h"

///////////////////////////////////////////////////////////////////
//
// Function Name:	Workstation
// Description:		Default constructor with no arguments.
//
///////////////////////////////////////////////////////////////////
Workstation::Workstation()
{
	active = false;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	Workstation
// Description:		Default destructor with no arguments.
//
///////////////////////////////////////////////////////////////////
Workstation::~Workstation()
{

}
