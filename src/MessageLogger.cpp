// ____________________________________________________________________________
//
//  General Information:
//
//  File Name:      MessageLogger.h
//  Author:         Timothy Hahn, Montana State University
//  Project:        RWASimulator
//
//  Description:    The file contains the implementation of the MessageLogger class,
//					used to record events to a file.
//
//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//  Revision History:
//
//  05/20/2009	v1.0	Initial Version.
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
	pthread_mutex_init(&LogMutex,NULL);
	pthread_mutex_init(&PrintMutex,NULL);
	pthread_mutex_init(&ResultsMutex,NULL);

	char buffer[100];
	sprintf(buffer,"output/EventLog-%s-%s-%s-%s-R%d.txt",topo,lambda,seed,k,runCount);

	eventLogger.open(buffer);
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
void MessageLogger::recordEvent(const std::string &e, bool print, unsigned short int ci)
{
	if(threadZero->getQualityParams().detailed_log == false && print == false)
	{
		return;
	}

	pthread_mutex_lock(&LogMutex);

	struct tm *current;
	time_t now;

	time(&now);
	current = localtime(&now);

	std::ostringstream message;

	char buff[25];

	strftime(buff, sizeof(buff), "%H:%I:%M", current);

	message << buff << " [] " << e << std::endl;

	eventLogger << message.str();

	pthread_mutex_unlock(&LogMutex);

	if(print == true)
	{
		pthread_mutex_lock(&PrintMutex);

		cout << message.str();

		pthread_mutex_unlock(&PrintMutex);
	}

}
