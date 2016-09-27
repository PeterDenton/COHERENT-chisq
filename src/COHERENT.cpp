#include <cmath>
#include <vector>
#include <cassert>
#include <fstream>

#include "COHERENT.h"
#include "State.h"
#include "Constants.h"
#include "Target.h"

#include <iostream>

double h_sigma = 0.1; // start with 10% from subsection A on page 4 of http://arxiv.org/pdf/hep-ex/0511042.pdf

double Qwalphasq(int alpha, Target t, State s)
{
	double gVp = 0.5 - 2 * ssqw;
	double gVn = -0.5;

	double Qwalphasq = 0;
	double epsaau = 0;
	double epsaad = 0;
	double epsa1u = 0;
	double epsa1d = 0;
	double epsa2u = 0;
	double epsa2d = 0;

	switch (alpha)
	{
		case 0:
			if (s.up)
			{
				epsaau = s.ee;
				epsa1u = s.em;
				epsa2u = s.et;
			}
			else
			{
				epsaad = s.ee;
				epsa1d = s.em;
				epsa2d = s.et;
			}
			break;
		case 1:
			if (s.up)
			{
				epsaau = s.mm;
				epsa1u = s.em;
				epsa2u = s.mt;
			}
			else
			{
				epsaad = s.mm;
				epsa1d = s.em;
				epsa2d = s.mt;
			}
			break;
		case 2:
			if (s.up)
			{
				epsaau = s.tt;
				epsa1u = s.et;
				epsa2u = s.mt;
			}
			else
			{
				epsaad = s.tt;
				epsa1d = s.et;
				epsa2d = s.mt;
			}
			break;
	}

	Qwalphasq += pow(t.Z * (gVp + 2 * epsaau + epsaad) + t.N * (gVn + epsaau + 2 * epsaad), 2);
	Qwalphasq += pow(t.Z * (2 * epsa1u + epsa1d) + t.N * (epsa1u + 2 * epsa1d), 2);
	Qwalphasq += pow(t.Z * (2 * epsa2u + epsa2d) + t.N * (epsa2u + 2 * epsa2d), 2);
	return 4 * Qwalphasq;
}

// returns an energy, a delta function really
double Enumu()
{
	return (pow(mpi, 2) - pow(mmu, 2)) / (2 * mpi);
}
// the normalized spectra
double pnumubar(double Enu)
{
	double x = Enu / mmu;
	if (x < 0 or x > 0.5)
		return 0;
	return (64 / mmu) * pow(x, 2) * (0.75 - x); // 64/mmu = normalization factor
}
// the normalized spectra
double pnue(double Enu)
{
	double x = Enu / mmu;
	if (x < 0 or x > 0.5)
		return 0;
	return (192 / mmu) * pow(x, 2) * (0.5 - x); // 192/mmu = normalization factor
}

// see eq 1 of hep-ex/0511042
// E is the nuclear recoil energy (MeV)
// M is the nuclear mass (MeV)
// k is the nu momentum (energy) (MeV)
// returns cm^2/MeV
double dsigmadE(double E, Target t, double k)
{
	double Qw = t.N - (1 - 4 * ssqw) * t.Z;
	double dsigmadE = (pow(GF, 2) / (2 * M_PI)) * (pow(Qw, 2) / 4) * pow(FF(2 * t.M * E, t), 2) * t.M * (2 - t.M * E / pow(k, 2)); // MeV / GeV^4
//std::cout << E << " " << k << " " << 2 - t.M * E / pow(k, 2) << std::endl;
	return std::max(0., 0.3894 * 1e-33 * dsigmadE); // cm^2 / MeV
}
// see eq 2 of hep-ex/0511042
// E is the nuclear recoil energy
// M is the nuclear mass
// flavor:
//	0 - mu
//	1 - mubar
//	2 - e
// Mtarget in kg
// flux in N/cm^2
double dNdE(double E, Target t, int flavor, double Mtarget, double flux)
{
	double Nt = Mtarget * 1e3 * NA / (t.M / amu2MeV); // number of targets, 1e3 to convert kg to g, amu2MeV converts M back from MeV to amu

	if (flavor == 0)
		return Nt * flux * dsigmadE(E, t, Enumu()); // delta function is an easy integral

	int Nbins = 1e2;
	double dNdE = 0;
	double kmax = mmu / 2;
	double kstep = kmax / Nbins;
	double phi; // normalized flux at a given energy
	for (double k = kstep / 2; k <= kmax; k += kstep) // neutrino momentum (energy) (MeV); start halfway in
	{
		phi = (flavor == 1) ? pnumubar(k) : pnue(k);
		dNdE += kstep * Nt * flux * phi * dsigmadE(E, t, k);
	} // k, neutrino momentum (energy)
	return dNdE;
}

// target recoil threshold in MeV, likely either 0.01 (10 keV) or 0.05 (50 keV)
// returns in cm^2
double g(int flavor, Target t, double Etr)
{
	double g = 0;
	double Emax = 0.3; // 300 keV seems to be enough
	int Nbins = 1e2;
	double Estep = (Emax - Etr) / Nbins;
	double kmax = mmu / 2;
	double kstep = kmax / Nbins;
	double phi; // normalized flux at a given energy
	for (double E = Etr; E < Emax; E+= Estep)
	{
		if (flavor == 0)
			g += Estep * (pow(GF, 2) / (8 * M_PI)) * t.M * pow(FF(2 * t.M * E, t), 2) * std::max(0., 2 - t.M * E / pow(Enumu(), 2));
		else
		{
			for (double k = kstep / 2; k <= kmax; k += kstep) // neutrino momentum (energy) (MeV); start halfway in
			{
				phi = (flavor == 1) ? pnumubar(k) : pnue(k);
				g += kstep * Estep * (pow(GF, 2) / (8 * M_PI)) * t.M * phi * pow(FF(2 * t.M * E, t), 2) * std::max(0., 2 - t.M * E / pow(k, 2)); // max to ensure positive definite, units: MeV^2/GeV^4
			} // k
		} // flavor = delayed
	} // E
	g *= 0.3894 / 1e6; // converts MeV^2/GeV^4 to mb
	g *= 1e-27; // converts mb to cm^2
	return g;
}

double g_pres[3][4][2];
void g_pre_init()
{
	double Etr;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			Etr = (j == 0) ? 0.01 : 0.05;

			g_pres[i][0][j] = g(i, Ne, Etr);
			g_pres[i][1][j] = g(i, Ar, Etr);
			g_pres[i][2][j] = g(i, Ge, Etr);
			g_pres[i][3][j] = g(i, Xe, Etr);
		} // Etr
	} // flavor
}
void g_pre_write()
{
	std::ofstream data("data/g_pres.txt");
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			for (int k = 0; k < 2; k++)
			{
				data << g_pres[i][j][k] << std::endl;
			} // k, Etr
		} // j, target
	} // i, flavor
	data.close();
}
void g_pre_read()
{
	std::ifstream data("data/g_pres.txt");
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			for (int k = 0; k < 2; k++)
			{
				data >> g_pres[i][j][k];
			} // k, Etr
		} // j, target
	} // i, flavor
	data.close();
}

double g_pre(int flavor, Target t, double Etr)
{
	assert(Etr == 0.01 or Etr == 0.05); // only pre calced for 10, 50 keV

	if (t.name == "Ne")
		return (Etr == 0.01 ? g_pres[flavor][0][0] : g_pres[flavor][0][1]);
	if (t.name == "Ar")
		return (Etr == 0.01 ? g_pres[flavor][1][0] : g_pres[flavor][1][1]);
	if (t.name == "Ge")
		return (Etr == 0.01 ? g_pres[flavor][2][0] : g_pres[flavor][2][1]);
	if (t.name == "Xe")
		return (Etr == 0.01 ? g_pres[flavor][3][0] : g_pres[flavor][3][1]);

	return 1e20;
}

// Mtargetyr in kg*yr
// Etr in MeV
double chisq(State s, std::vector<Target> ts, double Mtargetyr, double Etr, double sys)
{
	State sSM = {};
	sSM.up = s.up;

	double NpNSI, NpSM, NdNSI, NdSM, NmuNSI, NmuSM, NeNSI, NeSM, NmubarNSI, NmubarSM, contamination, chisq, Phi, Nt, galpha;
	chisq = 0;
	contamination = 0.138; // the amount of delayed events in the prompt window
	Phi = 1e7 * 60 * 60 * 24 * 365.25; // in #/cm^2/yr
	for (uint i = 0; i < ts.size(); i++)
	{
		Nt = Mtargetyr * 1e3 * NA / (ts[i].M / amu2MeV); // 1e3: kg -> g, need to convert M from MeV back to amu

		galpha = g_pre(2, ts[i], Etr); // electron
		NeNSI = Qwalphasq(0, ts[i], s) * Nt * Phi * galpha;
		NeSM = Qwalphasq(0, ts[i], sSM) * Nt * Phi * galpha;

		galpha = g_pre(0, ts[i], Etr); // muon
		NmuNSI = Qwalphasq(1, ts[i], s) * Nt * Phi * galpha;
		NmuSM = Qwalphasq(1, ts[i], sSM) * Nt * Phi * galpha;

		galpha = g_pre(1, ts[i], Etr); // antimuon
		NmubarNSI = Qwalphasq(1, ts[i], s) * Nt * Phi * galpha;
		NmubarSM = Qwalphasq(1, ts[i], sSM) * Nt * Phi * galpha;

		NpNSI = NmuNSI + contamination * (NeNSI + NmubarNSI);
		NpSM = NmuSM + contamination * (NeSM + NmubarSM);

		NdNSI = (1 - contamination) * (NeNSI + NmubarNSI);
		NdSM = (1 - contamination) * (NeSM + NmubarSM);

		if (NpSM == 0)
		{
			if (NpNSI != NpSM)
				chisq += 1e20;
		}
		else
			chisq += pow(NpNSI - NpSM, 2) / (NpSM + pow(sys * NpSM, 2));

		if (NdSM == 0)
		{
			if (NdNSI != NpSM)
				chisq += 1e20;
		}
		else
			chisq += pow(NdNSI - NdSM, 2) / (NdSM + pow(sys * NdSM, 2));
	} // i, targets
	return chisq;
}
