// ____________________________________________________________________________
//
//  General Information:
//
//  File Name:      Workstation.h
//  Author:         Timothy Hahn, PhD
//  Project:        raptor
//
//  Description:    The file contains the declaration of the Workstation class.
//					The purpose of the Workstations are to generate
//traffic for 					the optical network.
//
//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//  Revision History:
//
//  05/20/2009	v1.0	Initial Version.
//  04/14/2019  v2.0    Reworked version based upon cmake and octave
//
// ____________________________________________________________________________

#ifndef WORKSTATION_H
#define WORKSTATION_H

class Workstation {
 public:
  Workstation();
  ~Workstation();

  inline unsigned int getParentRouterIndex() { return parentRouterIndex; };
  inline void setParentRouterIndex(unsigned int p) { parentRouterIndex = p; };

  inline void setActive(bool a) { active = a; };
  inline bool getActive() { return active; };

 private:
  unsigned int parentRouterIndex;

  bool active;
};

#endif
