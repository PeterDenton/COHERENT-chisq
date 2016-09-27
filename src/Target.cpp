#include <fstream>
#include <string>

#include "Target.h"

#include <iostream>
#include <cmath>

// masses from: https://www.ncsu.edu/chemistry/msf/pdf/IsotopicMass_NaturalAbundance.pdf
Target Ne = {"Ne", 10, 10, 20, 19.992440 * amu2MeV};
Target Ar = {"Ar", 18, 22, 40, 39.962383 * amu2MeV};
Target Ge = {"Ge", 32, 44, 76, 75.921403 * amu2MeV};
Target Xe = {"Xe", 54, 78, 132, 131.904154 * amu2MeV};

// read in form factor data for a given target
void read_in(Target &t)
{
	if (not t.is_read_in)
	{
		std::string fname = "data/FF/";
		fname += t.name + std::to_string(t.A) + ".FF";
		std::ifstream dataf(fname);

		pair line;
		while (dataf >> line.x >> line.y) // Q in fm^-1, F(Q)
		{
			line.x *= 1e15; // Q in m^-1
			line.x *= 1.97327e-7; // Q in eV
			line.x *= 1e-6; // Q in MeV
			t.FF.push_back(line);
		}

		t.is_read_in = true;
	}
}

// interpolate the FF data from the tables
double FF(double Qsq, Target t)
{
	read_in(t);

	double Q = sqrt(Qsq);

	if (Q <= t.FF[0].x)
		return 1; // low Q => F ~= 1

	for (uint i = 1; i < t.FF.size(); i++)
	{
		if (t.FF[i].x > Q)
		{
			double m = (t.FF[i].y - t.FF[i - 1].y) / (t.FF[i].x - t.FF[i - 1].x);
			return m * (Q - t.FF[i].x) + t.FF[i].y;
		}
	} // i

	return 0; // if Q > max in data file, return 0
}
