// ____________________________________________________________________________
//
//  General Information:
//
//  File Name:      MessageLogger.h
//  Author:         Timothy Hahn, PhD
//  Project:        raptor
//
//  Description:    The file contains the implementation of the MessageLogger class,
//					used to record events to a file.
//
//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//  Revision History:
//
//  05/20/2009	v1.0	Initial Version.
//  04/14/2019  v2.0    Reworked version based upon cmake and octave
//
// ____________________________________________________________________________

#include "MessageLogger.h"
#include "Thread.h"

#include <cstdlib>
#include <ctime>
#include <sstream>

extern Thread* threadZero;
extern Thread** threads;

///////////////////////////////////////////////////////////////////
//
// Function Name:	MessageLogger
// Description:		Default constructor
//
///////////////////////////////////////////////////////////////////
MessageLogger::MessageLogger(const char* topo, const char* lambda, const char* seed, const char* k,int runCount)
{
	pthread_mutex_init(&LogMutex,nullptr);
	pthread_mutex_init(&PrintMutex,nullptr);
	pthread_mutex_init(&ResultsMutex,nullptr);

	std::ostringstream buffer;
	buffer << "output/EventLog-" << topo << "-" << lambda << "-" << seed << "-" << k << "-R" << runCount << ".txt";

	eventLogger.open(buffer.str());
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	~MessageLogger
// Description:		Default destructor
//
///////////////////////////////////////////////////////////////////
MessageLogger::~MessageLogger()
{
	pthread_mutex_destroy(&LogMutex);
	pthread_mutex_destroy(&PrintMutex);
	pthread_mutex_destroy(&ResultsMutex);

	eventLogger.close();
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	recordEvent
// Description:		Writes the event to the log file with a timestamp.
//					Also prints the file to the console if the print
//					value is set to true (that is the default)
//
///////////////////////////////////////////////////////////////////
void MessageLogger::recordEvent(const std::string &e, bool print, unsigned int ci)
{
	if(threadZero->getQualityParams().detailed_log == false && print == false)
	{
		return;
	}

	pthread_mutex_lock(&LogMutex);

	struct tm *current = nullptr;
	time_t now;

	time(&now);
#ifdef _MSC_VER
	current = new tm;
	localtime_s(current,&now);
#else
	current = localtime(&now);
#endif //_MSC_VER

	std::ostringstream message;

	char buff[25];

	strftime(buff, sizeof(buff), "%H:%M:%S", current);

	message << buff << " [] " << e << std::endl;

	eventLogger << message.str();

	pthread_mutex_unlock(&LogMutex);

	if(print == true)
	{
		pthread_mutex_lock(&PrintMutex);

		std::cout << message.str();

		pthread_mutex_unlock(&PrintMutex);
	}

	delete current;
}

///////////////////////////////////////////////////////////////////
//
// Function Name:	flushLog
// Description:		Flushes the log
//
///////////////////////////////////////////////////////////////////
void MessageLogger::flushLog(bool print)
{
	if (print == true)
	{
		pthread_mutex_lock(&PrintMutex);

		std::cout << std::flush;

		pthread_mutex_unlock(&PrintMutex);
	}
}
