#ifndef Target_H
#define Target_H

#include <string>
#include <vector>

#include "Constants.h"

struct pair
{
	double x, y;
}; // struct pair, used in FF

struct Target
{
	std::string name;
	int Z, N, A;
	double M; // MeV
	bool is_read_in = false; // whether or not the form factor has been read in
	std::vector<pair> FF; // form factor
}; // struct Target
extern Target Ne, Ar, Ge, Xe;

double FF(double Qsq, Target t);

#endif
