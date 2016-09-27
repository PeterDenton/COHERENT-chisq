#ifndef COHERENT_H
#define COHERENT_H

#include <vector>

#include "State.h"
#include "Target.h"

// Assumes that the measured central value is consistent with the SM
double Qwalphasq(int alpha, Target t, State s);
double chisq(State s, std::vector<Target> ts, double Mtargetyr, double Etr, double sys); // includes information about different flavors, Mtargetyr in kg*yr, sys is 0.1 for 10%

// the energy of the nu_mu
double Enumu();

// number of events as a function of nuclear recoil energy
double dNdE(double E, Target t, int flavor, double Mtarget, double flux);

extern double g_pres[3][4][2];
void g_pre_init();
void g_pre_write();
void g_pre_read();
double g_pre(int flavor, Target t, double Etr);

#endif
