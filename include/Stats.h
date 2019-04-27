// ____________________________________________________________________________
//
//  General Information:
//
//  File Name:      Stats.h
//  Author:         Timothy Hahn, PhD
//  Project:        raptor
//
//  Description:    The file contains datatypes that are important for
//  determining
//					the performance of the various
//algorithms.
//
//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//  Revision History:
//
//  05/20/2009	v1.0	Initial Version.
//  06/02/2009	v1.02	Minor optimizations and bug fixes.
//  04/14/2019  v2.0    Reworked version based upon cmake and octave
//
// ____________________________________________________________________________

#ifndef STATS_H
#define STATS_H

#include <cstdio>

#include "Edge.h"

struct GlobalStats {
  size_t ConnectionRequests;
  size_t ConnectionSuccesses;
  size_t CollisionFailures;
  size_t QualityFailures;
  size_t NoPathFailures;
  size_t DroppedFailures;
  size_t ProbeSentCount;
  size_t totalHopCount;
  size_t totalSpanCount;
  double aseNoiseTotal;
  double xpmNoiseTotal;
  double fwmNoiseTotal;
  double totalSetupDelay;
  double raRunTime;
};

struct EdgeStats {
  size_t droppedConnections;
  double minInitalQFactor;
  double minAverageQFactor;
  double minPercentQFactor;
  double maxInitalQFactor;
  double maxAverageQFactor;
  double maxPercentQFactor;
  double totalInitalQFactor;
  double totalAverageQFactor;
  double totalPercentQFactor;
  double totalTime;
  double count;
};

enum FailureTypes {
  COLLISION_FAILURE = -1,
  QUALITY_FAILURE = -2,
  NO_PATH_FAILURE = -3
};

#endif
