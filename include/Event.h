// ____________________________________________________________________________
//
//  General Information:
//
//  File Name:      Event.h
//  Author:         Timothy Hahn, PhD
//  Project:        raptor
//
//  Description:    The file contains the declaration of some data types
//					required by the EventQueue.
//
//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//  Revision History:
//
//  05/20/2009	v1.0	Initial Version.
//  06/02/2009	v1.02	Minor optimizations and bug fixes.
//  04/14/2019  v2.0    Reworked version based upon cmake and octave
//
// ____________________________________________________________________________

#ifndef EVENT_H
#define EVENT_H

#include "KShortestPaths.h"

#include "Edge.h"

#include <vector>

enum EventType {
	ACTIVATE_WORKSTATIONS,
	DEACTIVATE_WORKSTATIONS,
	UPDATE_USAGE,
	UPDATE_GUI,
	CONNECTION_REQUEST,
	CREATE_CONNECTION_PROBE,
	CREATE_CONNECTION_CONFIRMATION,
	COLLISION_NOTIFICATION,
	DESTROY_CONNECTION_PROBE,
	NUMBER_OF_EVENTS
};

struct Event {
	EventType e_type;
	double e_time;
	void* e_data;
};

struct ConnectionRequestEvent
{
	unsigned int sourceRouterIndex;
	unsigned int destinationRouterIndex;
	int wavelength;
	double connectionDuration;
	double requestBeginTime;
	unsigned int session;
	unsigned int sequence;
	unsigned int max_sequence;
	bool qualityFail;
};

struct CreateConnectionProbeEvent
{
	unsigned int sourceRouterIndex;
	unsigned int destinationRouterIndex;
	double connectionDuration;
	double requestBeginTime;
	double decisionTime;
	Edge **connectionPath;
	unsigned int connectionLength;
	unsigned int numberOfHops;
	int wavelength;
	unsigned int session;
	unsigned int sequence;
	unsigned int max_sequence;
	kShortestPathReturn *kPaths;
	CreateConnectionProbeEvent **probes;
	bool atDestination;
	bool qualityFail;
};

struct CollisionNotificationEvent
{
	unsigned int sourceRouterIndex;
	unsigned int destinationRouterIndex;
	Edge **connectionPath;
	unsigned int connectionLength;
	unsigned int numberOfHops;
	unsigned int session;
	unsigned int sequence;
	unsigned int max_sequence;
	int wavelength;
	CreateConnectionProbeEvent **probes;
	bool finalFailure;
};

struct CreateConnectionConfirmationEvent
{
	unsigned int sourceRouterIndex;
	unsigned int destinationRouterIndex;
	double connectionDuration;
	double requestBeginTime;
	Edge **connectionPath;
	unsigned int connectionLength;
	unsigned int numberOfHops;
	unsigned int max_sequence;
	int wavelength;
	int originalWavelength;
	unsigned int session;
	unsigned int sequence;
	kShortestPathReturn *kPaths;
	CreateConnectionProbeEvent **probes;
	bool finalFailure;
};

struct DestroyConnectionProbeEvent
{
	Edge **connectionPath;
	unsigned int connectionLength;
	unsigned int numberOfHops;
	unsigned int session;
	unsigned int sequence;
	int wavelength;
	CreateConnectionProbeEvent **probes;
};

#endif
