// ____________________________________________________________________________
//
//  General Information:
//
//  File Name:      ErrorCodes.h
//  Author:         Timothy Hahn, PhD
//  Project:        raptor
//
//  Description:    The file defines the various error codes returned by the
//					program at a critical error.
//
//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//  Revision History:
//
//  05/20/2009	v1.0	Initial Version.
//  04/14/2019  v2.0    Reworked version based upon cmake and octave
//
// ____________________________________________________________________________

#ifndef ERROR_CODES_H
#define ERROR_CODES_H

enum ErrorCodes {
  ERROR_INVALID_PARAMETERS = -1,
  ERROR_TOO_MANY_EDGES = -2,
  ERROR_CHOOSE_WAVELENGTH_1 = -3,
  ERROR_CHOOSE_WAVELENGTH_2 = -4,
  ERROR_INVALID_EDGE = -5,
  ERROR_THREAD_INIT = -6,
  ERROR_THREAD_EVENT_TYPE = -7,
  ERROR_EDGE_IS_USED = -8,
  ERROR_EDGE_IS_FREE = -9,
  ERROR_INVALID_CONFIRMATION = -10,
  ERROR_QUALITY_INPUT = -11,
  ERROR_TOPOLOGY_INPUT_ROUTERS = -12,
  ERROR_TOPOLOGY_INPUT_EDGES = -13,
  ERROR_WORKSTATION_INPUT_QUANTITY = -14,
  ERROR_WORKSTATION_INPUT_PARENT = -15,
  ERROR_ALGORITHM_INPUT = -16,
  ERROR_RECORD_EVENT = -17,
  ERROR_PROBES_RECEIVED = -18,
  ERROR_USER_CLOSED = -19,
  ERROR_QFACTOR_MONTIOR = -20,
  ERROR_WAVELENGTH_ALGORITHM_IA = -21,
  ERROR_PRIORITY_QUEUE = -22,
  ERROR_WORKSTATION_FILE = -23,
  ERROR_QUALITY_FILE = -24,
  ERROR_TOPOLOGY_FILE = -25,
  ERROR_NO_FLUSH = -26,
  ERROR_GUI = -27,
  ERROR_THREAD_CREATION = -28,
  ERROR_OCTAVE = -29
};

#endif
