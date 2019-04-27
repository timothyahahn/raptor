// ____________________________________________________________________________
//
//  General Information:
//
//  File Name:      Edge.h
//  Author:         Timothy Hahn, PhD
//  Project:        raptor
//
//  Description:    The file contains the declaration of the Edge class.
//					The purpose of the Edge class is to connect the
//routers 					together.
//
//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//  Revision History:
//
//  05/20/2009	v1.0	Initial Version.
//  06/02/2009	v1.02	Minor optimizations and bug fixes.
//  04/14/2019  v2.0    Reworked version based upon cmake and octave
//
// ____________________________________________________________________________

#ifndef EDGE_H
#define EDGE_H

#include <list>
#include <vector>

#include "Stats.h"

#ifndef NO_ALLEGRO
#include "allegro5/allegro.h"
#endif

enum EdgeStatus { EDGE_FREE, EDGE_USED };

class Edge {
 public:
  Edge();
  Edge(size_t src, size_t dest, size_t spans);

  ~Edge();

  inline size_t getSourceIndex() const { return sourceIndex; }
  inline size_t getDestinationIndex() const { return destinationIndex; }
  inline size_t getNumberOfSpans() const { return numberOfSpans; }
  inline EdgeStatus getStatus(size_t w) const { return status[w]; }
  inline size_t getActiveSession(long long int w) const { return activeSession[w]; }
  inline double getAlgorithmUsage() const { return algorithmUsage; }
  inline double getQMDegredation() const { return QMDegredation; }
  inline double getPheremone() const { return pheremone; }
  inline EdgeStats* getEdgeStats() const { return stats; }

  inline void setUsed(size_t session, size_t w) {
    status[w] = EDGE_USED;
    activeSession[w] = session;
  };
  inline void setFree(long long int w) {
    status[w] = EDGE_FREE;
    activeSession[w] = -1;
    degredation[w] = 0.0;
  };

  void updateUsage();

  void updateQMDegredation(size_t ci, long long int wavelength);
  void updateQFactorStats(size_t ci, long long int wavelength);

  inline void resetAlgorithmUsage() { algorithmUsage = 0.0; };
  void resetQMDegredation();

  void resetPheremone(size_t ci, size_t spans);

  void evaporatePheremone(size_t ci);
  void addPheremone(size_t hops, size_t ci);

#ifndef NO_ALLEGRO
  void updateGUI();
  void refreshbmps(bool useThread);
  void initializetopobmps();
  inline int getMaxActualUsage() { return max_actual_usage; };
  inline int getUsageNums() { return static_cast<int>(usageList.size()); };
  inline void addUsage(int u) { usageList.push_back(u); }
  inline void setMaxUsage(int u) { max_actual_usage = u; }
  void scaleEdgesTo(int spns, int px);
  void paintSpans();
  void paintUsage(int p);
  void paintUsageHistory(int x1, int y1, int x2, int y2, int t);
  void calculateAverageUsage();
  void saveData(char* file);
#endif
  void resetEdgeStats();

  inline void insertEstablishedConnection(void* ec_void) {
    establishedConnections.push_back(ec_void);
  };
  void removeEstablishedConnection(void* dcpe_void);

 private:
  size_t sourceIndex;
  size_t destinationIndex;
  size_t numberOfSpans;

  long long int* activeSession;
  EdgeStatus* status;

  double algorithmUsage;
  int actualUsage;

  double QMDegredation;

  double* degredation;

  std::list<void*> establishedConnections;

#ifndef NO_ALLEGRO
  int r1x, r1y, r2x, r2y;  // coordinates of routers
  int r3x, r3y, r4x, r4y;  // points to draw to for edge width (these change)
  int putX, putY;
  int thdIndx;
  double invunitX, invunitY, currX, currY, oldX, oldY, average_usage,
      painted_usage;  // inverse unit vectors
  int painted_percent;
  int max_actual_usage;
  std::vector<int> usageList;
  ALLEGRO_BITMAP* edgeBmps[14];
  ALLEGRO_BITMAP* edgeBmp;
#endif
  EdgeStats* stats;

  double pheremone;
};

#endif
