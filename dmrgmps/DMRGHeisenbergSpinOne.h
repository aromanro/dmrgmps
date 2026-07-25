#pragma once
#include "MPSDMRGAlgorithm.h"
#include "SiteOperator.h"

namespace DMRG {

	namespace Heisenberg {

		class DMRGHeisenbergSpinOne :
			public MPS::MPSDMRGAlgorithm
		{
		public:
			DMRGHeisenbergSpinOne(double Jz = 1., double Jxy = 1., unsigned int maxstates = 10, unsigned int nrStates = 0, unsigned int method = 0)
				: MPS::MPSDMRGAlgorithm(3, Operators::SzOne(3).matrix, Operators::SplusOne(3).matrix, Jz, Jxy, maxstates, nrStates, method)
			{}
		};

	}
}

