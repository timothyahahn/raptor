// ____________________________________________________________________________
//
//  General Information:
//
//  File Name:      MessageLogger.h
//  Author:         Timothy Hahn, PhD
//  Project:        raptor
//
//  Description:    The file contains the declaration of the MessageLogger class,
//					used to record events to a file.
//
//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//  Revision History:
//
//  05/20/2009	v1.0	Initial Version.
//  04/14/2019  v2.0    Reworked version based upon cmake and octave
//
// ____________________________________________________________________________

#ifndef MESSAGE_LOGGER_H
#define MESSAGE_LOGGER_H

#define HAVE_STRUCT_TIMESPEC
#include "pthread.h"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>

class MessageLogger
{
	public:
		MessageLogger(const std::string&, const std::string&, const std::string&, const std::string&);
		~MessageLogger();

		void recordEvent(const std::string &e, bool print, unsigned int ci);

		void flushLog(bool print);

		inline void LockResultsMutex()
			{ pthread_mutex_lock(&ResultsMutex); };

		inline void UnlockResultsMutex()
			{ pthread_mutex_unlock(&ResultsMutex); };

	private:
		std::ofstream eventLogger;

		pthread_mutex_t LogMutex;
		pthread_mutex_t PrintMutex;
		pthread_mutex_t ResultsMutex;
};

#endif
