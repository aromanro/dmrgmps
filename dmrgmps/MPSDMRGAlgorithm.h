#pragma once

// Variational single-site MPS-DMRG ground state solver.
//
// This replaces the old block / reduced-density-matrix DMRG. The ground state is
// searched directly as a Matrix Product State: the chain is swept left/right and,
// at each site, the effective single-site Hamiltonian (built from the MPO and the
// left/right environments) is diagonalized with Lanczos. The orthogonality center
// is moved with a truncated SVD (bond dimension limited by 'maxStates').
//
// The public interface matches what the GUI / DMRGThread expect:
//   - constructor style used by the concrete Heisenberg classes,
//   - nrStates, EnergyGap, results, truncationError,
//   - CalculateFinite / CalculateInfinite.

#include "MPS.h"
#include "EffectiveHamiltonian.h"

#include <list>
#include <vector>

namespace DMRG {

	namespace MPS {

		// How the orthogonality center is optimized while sweeping.
		enum class SweepMode
		{
			SingleSite = 0,             // one-site update (fastest, can get stuck in local minima)
			SingleSiteSubspaceExpansion = 1, // one-site update + density-matrix (White) perturbation
			TwoSite = 2                 // two-site update (slower, allows bond growth, more robust)
		};

		class MPSDMRGAlgorithm
		{
		public:
			MPSDMRGAlgorithm(int physicalDim, const Tensor2& Sz, const Tensor2& Splus, double Jz = 1., double Jxy = 1., unsigned int maxstates = 10, unsigned int nrStates = 0, unsigned int method = 0);
			virtual ~MPSDMRGAlgorithm() = default;

			// results reported to the GUI
			double truncationError;
			std::list<double> results; // per-bond energy contribution <h_i,i+1> of the final ground state
			unsigned int nrStates;     // number of excited states targeted (for the gap)
			double EnergyGap;
			double groundEnergy;

			// left/right MPO environments per site
			std::vector<Tensor3> Lenv;
			std::vector<Tensor3> Renv;
			// overlap environments per reference (only used for deflation)
			std::vector<std::vector<Tensor2>> OLenv;
			std::vector<std::vector<Tensor2>> ORenv;

			// energies of the excited states found (size nrStates when nrStates > 0)
			std::vector<double> excitedEnergies;

			// --- solver options ------------------------------------------------
			// which sweep update to use
			SweepMode sweepMode = SweepMode::SingleSite;
			// density-matrix perturbation amplitude for SingleSiteSubspaceExpansion
			// (mixed in on every sweep except the last, where it is set to zero)
			double perturbationFactor = 1E-2;
			// base RNG seed for the (reproducible) random initial state; excited
			// states are seeded with rngSeed + stateIndex so they do not start out
			// identical to the ground state
			unsigned int rngSeed = 42;

			double CalculateFinite(int chainLength, int numSweeps);

		protected:
			int m_d;
			Tensor2 m_Sz;
			Tensor2 m_Splus;
			double m_Jz;
			double m_Jxy;
			unsigned int maxStates;
			double svdThreshold;

			HeisenbergMPO mpo;
			MatrixProductState mps;

			void SweepRightTwoSites(MatrixProductState& state, const std::vector<const MatrixProductState*>& references, double& energy, double& maxTruncation, bool deflate, double penaltyWeight);
			void SweepLeftTwoSites(MatrixProductState& state, const std::vector<const MatrixProductState*>& references, double& energy, double& maxTruncation, bool deflate, double penaltyWeight);
			void SweepRightOneSite(MatrixProductState& state, const std::vector<const MatrixProductState*>& references, double& energy, double& maxTruncation, bool deflate, double penaltyWeight, double alpha);
			void SweepLeftOneSite(MatrixProductState& state, const std::vector<const MatrixProductState*>& references, double& energy, double& maxTruncation, bool deflate, double penaltyWeight, double alpha);

			Tensor4 SolvePair(int i, MatrixProductState& state, const std::vector<const MatrixProductState*>& references, double& energy, bool deflate, double penaltyWeight);
			void SolveSite(int i, MatrixProductState& state, const std::vector<const MatrixProductState*>& references, double& energy, bool deflate, double penaltyWeight);

			Eigen::MatrixXd BuildSitePenalty(int i, int Dl, int Dr, const std::vector<const MatrixProductState*>& references, bool deflate) const;

			// Runs the sweeps to relax 'state' towards the (possibly deflated) ground
			// state. Each state in 'references' contributes a penalty
			// weight * |ref><ref| so the search converges to a state orthogonal to all
			// of them (used to target excited states). The references only enter via
			// gauge-independent overlaps, so their canonical form does not matter.
			double GroundStateSearch(MatrixProductState& state, int numSweeps, const std::vector<const MatrixProductState*>& references, double penaltyWeight = 0.);

			double EnergyExpectation(const MatrixProductState& state);
			void ComputeMeasurements(MatrixProductState& state);
		};

	}

}
