#include "stdafx.h"

#define _DMRGThread

#include "DMRGThread.h"


template<class Algorithm> DMRGThread<Algorithm>::DMRGThread(int sites, double Jz, double Jxy, int sweeps, int states, int nrExcitedStates, unsigned int method)
	: dmrg(Jz, Jxy, states, nrExcitedStates, method), m_Sites(sites), m_Sweeps(sweeps)
{
}

template<class Algorithm> void DMRGThread<Algorithm>::Calculate()
{
	result = dmrg.CalculateFinite(m_Sites, m_Sweeps);
	gapResult = dmrg.EnergyGap;

	terminated = true;
}