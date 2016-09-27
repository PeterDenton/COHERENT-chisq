#include <vector>
#include <fstream>

#include "Target.h"
#include "COHERENT.h"

int main()
{
//	g_pre_init(); // integrate the g functions (takes some time)
//	g_pre_write(); // write the g functions to file
	g_pre_read(); // read in the pre integrated g functions

	double Etr = 0.01; // MeV, works for either 0.01 or 0.05. 0.01 is optimistic, 0.05 is likely
	double sys = 0; // fractional, so 0, 0.05, 0.1 are reasonable values
	double Mtargetyr = 100; // the mass of the target in kg times the number of years (100 kg*yr is reasonable)

	std::vector<Target> ts; // the list of targets
	ts.push_back(Ar); // Argon
	ts.push_back(Xe); // Xenon

	State s = {}; // the NSI object, defaults to the SM
	s.up = true; // since we are only doing one quark at a time, we must specify which is 'on.' True => up, false => down.

	std::ofstream data("data/ee.txt"); // file to write to

	for (double ee = -0.1; ee <= 0.6; ee += 0.001)
	{
		s.ee = ee; // set the eps_ee parameter
		data << ee << " " << chisq(s, ts, Mtargetyr, Etr, sys) << std::endl; // write the chisq to file
	} // ee

	data.close();

	return 0;
}
