#pragma once
#include "MPSDMRGAlgorithm.h"
#include "SiteOperator.h"

namespace DMRG {

	namespace Heisenberg {

		class DMRGHeisenbergSpinOneHalf :
			public MPS::MPSDMRGAlgorithm
		{
		public:
			DMRGHeisenbergSpinOneHalf(double Jz = 1., double Jxy = 1., unsigned int maxstates = 18, unsigned int nrStates = 0, unsigned int method = 0)
				: MPS::MPSDMRGAlgorithm(2, Operators::SzOneHalf(2).matrix, Operators::SplusOneHalf(2).matrix, Jz, Jxy, maxstates, nrStates, method)
			{}
		};
	}
}