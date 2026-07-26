#include "MPSDMRGAlgorithm.h"

#include <algorithm>
#include <limits>
#include <cmath>
#include <cassert>

#include <iostream>


namespace DMRG {

	namespace MPS {

		namespace
		{
			constexpr double LANCZOS_EIGENVALUE_TOLERANCE = 1E-8;
			constexpr double LANCZOS_RESIDUAL_TOLERANCE = 1E-10;

			// overlap environment updates between a 'current' state and a fixed 'reference' state

			inline Tensor2 UpdateLeftOverlap(const Tensor2& OL, const Tensor3& cur, const Tensor3& ref)
			{
				// cur (i, j, k), OL (k, l), ref (n, j, l)
				const Eigen::array<Eigen::IndexPair<int>, 1> c1{ Eigen::IndexPair<int>(0, 0) };
				// first contraction: cur (i, j, k) OL (k, l) -> (i, j, l)

				const Eigen::array<Eigen::IndexPair<int>, 2> c2{ Eigen::IndexPair<int>(0, 1), Eigen::IndexPair<int>(2, 0) };
				// second contraction: (i, j, l) ref (n, j, l) -> (i, n)

				Tensor2 result = cur.contract(OL, c1).contract(ref, c2);

				assert (result.dimension(0) == cur.dimension(2) && result.dimension(1) == ref.dimension(2));

				return result;
			}

			inline Tensor2 UpdateRightOverlap(const Tensor2& OR, const Tensor3& cur, const Tensor3& ref)
			{
				// cur (i, j, k), OR (i, l), ref (l, j, n)
				const Eigen::array<Eigen::IndexPair<int>, 1> c1{ Eigen::IndexPair<int>(2, 0) };
				// first contraction: cur (i, j, k) OR (i, l) -> (j, k, l)

				const Eigen::array<Eigen::IndexPair<int>, 2> c2{ Eigen::IndexPair<int>(1, 1), Eigen::IndexPair<int>(2, 2) };
				// second contraction: (j, k, l) ref (l, j, n) -> (k, n)

				Tensor2 result = cur.contract(OR, c1).contract(ref, c2);

				assert (result.dimension(0) == cur.dimension(0) && result.dimension(1) == ref.dimension(0));

				return result;
			}

			// Full overlap <a|b> of two (unnormalized) MPS by folding the transfer
			// matrices from left to right. Both are assumed to have singleton left and
			// right boundary bonds.
			inline double Overlap(const MatrixProductState& a, const MatrixProductState& b)
			{
				Tensor2 O(1, 1);
				O(0, 0) = 1.;
				for (int i = 0; i < a.N; ++i)
					O = UpdateLeftOverlap(O, a.GetSiteTensor(i), b.GetSiteTensor(i));

				assert(O.dimension(0) == 1 && O.dimension(1) == 1);

				return O(0, 0);
			}

			// Density-matrix (White / Hubig) subspace-expansion mixing blocks.
			//
			// Left-to-right: build P(l', s', b, r) = alpha * sum_{a,l,s} L(l,a,l')
			//   W(a,b,s',s) A(l,s,r), reshaped into a (Dl*d) x (Dw*Dr) block of extra
			//   columns appended to the left-normalization matrix.
			Eigen::MatrixXd LeftMixing(const HeisenbergMPO& mpo, const Tensor3& L, const Tensor3& A, double alpha)
			{
				const int Dl = static_cast<int>(A.dimension(0));
				const int dd = static_cast<int>(A.dimension(1));
				const int Dr = static_cast<int>(A.dimension(2));
				const int Dw = HeisenbergMPO::Dw;

				const Eigen::array<Eigen::IndexPair<int>, 1> c1{ Eigen::IndexPair<int>(0, 0) };
				const Eigen::Tensor<double, 4> T1 = L.contract(A, c1); // (a, l', s, r)

				const Eigen::array<Eigen::IndexPair<int>, 2> c2{ Eigen::IndexPair<int>(0, 0), Eigen::IndexPair<int>(2, 3) };
				const Eigen::Tensor<double, 4> P = T1.contract(mpo.W, c2); // (l', r, b, s')

				Eigen::MatrixXd M(static_cast<Eigen::Index>(Dl) * dd, static_cast<Eigen::Index>(Dw) * Dr);
				for (int lp = 0; lp < Dl; ++lp)
					for (int r = 0; r < Dr; ++r)
						for (int b = 0; b < Dw; ++b)
							for (int sp = 0; sp < dd; ++sp)
								M(static_cast<Eigen::Index>(sp) * Dl + lp, static_cast<Eigen::Index>(b) * Dr + r) = alpha * P(lp, r, b, sp);

				return M;
			}

			// Right-to-left: build P(a, l, s', r') = alpha * sum_{b,r,s} R(r,b,r')
			//   W(a,b,s',s) A(l,s,r), reshaped into a (Dw*Dl) x (d*Dr) block of extra
			//   rows appended to the right-normalization matrix.
			Eigen::MatrixXd RightMixing(const HeisenbergMPO& mpo, const Tensor3& R, const Tensor3& A, double alpha)
			{
				const int Dl = static_cast<int>(A.dimension(0));
				const int dd = static_cast<int>(A.dimension(1));
				const int Dr = static_cast<int>(A.dimension(2));
				const int Dw = HeisenbergMPO::Dw;

				const Eigen::array<Eigen::IndexPair<int>, 1> c1{ Eigen::IndexPair<int>(0, 2) };
				const Eigen::Tensor<double, 4> T1 = R.contract(A, c1); // (b, r', l, s)

				const Eigen::array<Eigen::IndexPair<int>, 2> c2{ Eigen::IndexPair<int>(0, 1), Eigen::IndexPair<int>(3, 3) };
				const Eigen::Tensor<double, 4> P = T1.contract(mpo.W, c2); // (r', l, a, s')

				Eigen::MatrixXd M(static_cast<Eigen::Index>(Dw) * Dl, static_cast<Eigen::Index>(dd) * Dr);
				for (int rp = 0; rp < Dr; ++rp)
					for (int l = 0; l < Dl; ++l)
						for (int a = 0; a < Dw; ++a)
							for (int sp = 0; sp < dd; ++sp)
								M(static_cast<Eigen::Index>(a) * Dl + l, static_cast<Eigen::Index>(sp) * Dr + rp) = alpha * P(rp, l, a, sp);

				return M;
			}

			// Lanczos ground state of any object exposing Size() and Apply(const VectorXd&).
			// Uses full reorthogonalization; the Krylov basis is grown lazily in blocks so
			// that quick convergence does not pay for a full n x maxIter allocation.
			template <typename Ham>
			double LanczosGroundState(const Ham& H, Eigen::VectorXd& groundState)
			{
				const Eigen::Index n = H.Size();

				if (n == 1)
				{
					const Eigen::VectorXd v = Eigen::VectorXd::Ones(1);
					const double e = H.Apply(v)(0);
					groundState = v;
					return e;
				}

				const int maxIter = static_cast<int>(std::min<Eigen::Index>(n, 150));
				const int growBlock = static_cast<int>(std::min<Eigen::Index>(maxIter, 16));

				Eigen::VectorXd q;
				if (groundState.size() == n && groundState.norm() > 1E-12)
					q = groundState.normalized();
				else
				{
					q = Eigen::VectorXd::Random(n);
					q.normalize();
				}

				Eigen::MatrixXd basis(n, std::min(growBlock, maxIter));
				Eigen::VectorXd alpha = Eigen::VectorXd::Zero(maxIter);
				Eigen::VectorXd beta = Eigen::VectorXd::Zero(maxIter);
				Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver;

				basis.col(0) = q;
				double groundEnergy = std::numeric_limits<double>::infinity();
				int krylovDim = 1;

				for (int j = 0; j < maxIter; ++j)
				{
					q = basis.col(j);
					Eigen::VectorXd z = H.Apply(q);
					if (j > 0)
						z -= beta(j - 1) * basis.col(j - 1);

					alpha(j) = q.dot(z);
					z -= alpha(j) * q;

					for (int k = 0; k <= j; ++k)
						z -= basis.col(k).dot(z) * basis.col(k);

					krylovDim = j + 1;
					if (krylovDim > 3)
					{
						solver.computeFromTridiagonal(alpha.head(krylovDim), beta.head(krylovDim - 1));
						const double eigenval = solver.eigenvalues()(0);
						if (std::abs(groundEnergy - eigenval) < LANCZOS_EIGENVALUE_TOLERANCE * std::max(1., std::abs(eigenval)))
						{
							groundEnergy = eigenval;
							break;
						}

						groundEnergy = eigenval;
					}

					if (j == maxIter - 1)
						break;

					beta(j) = z.norm();
					if (beta(j) < LANCZOS_RESIDUAL_TOLERANCE)
						break;

					// grow the basis lazily, making room for the next Krylov vector
					if (j + 1 >= basis.cols())
						basis.conservativeResize(Eigen::NoChange, std::min(maxIter, static_cast<int>(basis.cols()) + growBlock));

					basis.col(j + 1) = z / beta(j);
				}

				if (krylovDim == 1)
				{
					groundState = basis.col(0);
					return alpha(0);
				}

				solver.computeFromTridiagonal(alpha.head(krylovDim), beta.head(krylovDim - 1));
				groundEnergy = solver.eigenvalues()(0);
				groundState = basis.leftCols(krylovDim) * solver.eigenvectors().col(0);

				const double nrm = groundState.norm();
				if (nrm > 0.) groundState /= nrm;

				return groundEnergy;
			}
		}


		MPSDMRGAlgorithm::MPSDMRGAlgorithm(int physicalDim, const Tensor2& Sz, const Tensor2& Splus, double Jz, double Jxy, unsigned int maxstates, unsigned int nrStates, unsigned int method)
			: truncationError(0), nrStates(nrStates), EnergyGap(0), groundEnergy(0), sweepMode(static_cast<SweepMode>(method)),
			m_d(physicalDim), m_Sz(Sz), m_Splus(Splus), m_Jz(Jz), m_Jxy(Jxy), maxStates(maxstates), svdThreshold(1E-10)
		{
			mpo.Build(m_d, m_Sz, m_Splus, m_Jz, m_Jxy);
		}

		double MPSDMRGAlgorithm::CalculateFinite(int chainLength, int numSweeps)
		{
			if (chainLength <= 0) return std::numeric_limits<double>::infinity();
			if (numSweeps < 1) numSweeps = 1;

			results.clear();
			excitedEnergies.clear();
			truncationError = 0;
			EnergyGap = 0;

			const int chi = static_cast<int>(maxStates);

			mps.Init(chainLength, m_d, chi, rngSeed);
			mps.RightCanonicalize(chi, svdThreshold);

			groundEnergy = GroundStateSearch(mps, numSweeps, {}, 0.);

			// per-site measurements of the ground state
			ComputeMeasurements(mps);

			// optional excited states: each new state is deflated against the ground
			// state and against all previously found excited states. Overlaps are
			// verified after the search; if a state is not sufficiently orthogonal the
			// penalty weight is increased and the search is retried.
			if (nrStates > 0)
			{
				const double baseWeight = 10. * (std::abs(groundEnergy) + 1.);
				const double overlapTol = 1E-3;

				// keep the excited states alive so later ones can deflate against them
				std::vector<MatrixProductState> found;
				found.reserve(nrStates);

				for (unsigned int k = 0; k < nrStates; ++k)
				{
					std::vector<const MatrixProductState*> references;
					references.reserve(1 + found.size());
					references.push_back(&mps);
					for (const auto& s : found)
						references.push_back(&s);

					MatrixProductState excited;
					double excitedEnergy = std::numeric_limits<double>::quiet_NaN();

					double penaltyWeight = baseWeight;
					for (int attempt = 0; attempt < 3; ++attempt)
					{
						excited.Init(chainLength, m_d, chi, rngSeed + k + 11 + static_cast<unsigned int>(attempt) * 107u);
						excited.RightCanonicalize(chi, svdThreshold);

						GroundStateSearch(excited, numSweeps, references, penaltyWeight);
						excitedEnergy = EnergyExpectation(excited);

						double maxOverlap = 0.;
						for (const auto* ref : references)
							maxOverlap = std::max(maxOverlap, std::abs(Overlap(excited, *ref)));

						if (maxOverlap < overlapTol)
							break;

						penaltyWeight *= 10.;
					}

					excitedEnergies.push_back(excitedEnergy);
					found.push_back(std::move(excited));
				}

				const double lastExcited = excitedEnergies.back();
				//const double previousExcited = (excitedEnergies.size() > 1) ? excitedEnergies[excitedEnergies.size() - 2] : groundEnergy;
				const double gap = lastExcited - groundEnergy;
				if (std::isfinite(lastExcited) && std::isfinite(gap) && gap >= -LANCZOS_EIGENVALUE_TOLERANCE * std::max(1., std::abs(groundEnergy)))
					EnergyGap = std::max(0., gap);
				else
					EnergyGap = std::numeric_limits<double>::quiet_NaN();
			}

			return groundEnergy;
		}

		void MPSDMRGAlgorithm::InitREnv(const MatrixProductState& state)
		{
			const int N = state.N;

			Renv.clear();
			Renv.resize(N);

			Renv[N - 1] = mpo.RightBoundary();
			for (int i = N - 1; i > 0; --i)
				Renv[i - 1] = mpo.UpdateRight(Renv[i], state.GetSiteTensor(i));
		}

		double MPSDMRGAlgorithm::GroundStateSearch(MatrixProductState& state, int numSweeps, const std::vector<const MatrixProductState*>& references, double penaltyWeight)
		{
			const int N = state.N;
			const bool deflate = !references.empty() && (penaltyWeight != 0.);
			const int R = deflate ? static_cast<int>(references.size()) : 0;
			const bool twoSite = (sweepMode == SweepMode::TwoSite) && (N >= 2);
			const bool useMixing = (sweepMode == SweepMode::SingleSiteSubspaceExpansion) && (perturbationFactor > 0.) && !twoSite;

			// left/right MPO environments per site
			Lenv.clear();
			Lenv.resize(N);
			Lenv[0] = mpo.LeftBoundary();

			// state is right canonical (center at site 0): build all right environments
			InitREnv(state);

			// overlap environments per reference (only used for deflation)
			OLenv.clear();
			OLenv.resize(R, std::vector<Tensor2>(N));
			ORenv.clear();
			ORenv.resize(R, std::vector<Tensor2>(N));

			for (int r = 0; r < R; ++r)
			{
				OLenv[r][0] = Tensor2(1, 1);
				OLenv[r][0](0, 0) = 1.;
				ORenv[r][N - 1] = Tensor2(1, 1);
				ORenv[r][N - 1](0, 0) = 1.;
				for (int i = N - 1; i > 0; --i)
					ORenv[r][i - 1] = UpdateRightOverlap(ORenv[r][i], state.GetSiteTensor(i), references[r]->GetSiteTensor(i));
			}

			double energy = std::numeric_limits<double>::infinity();
			double maxTruncation = 0.;

			if (twoSite)
			{
				for (int sweep = 0; sweep < numSweeps; ++sweep)
				{
					// left to right
					SweepRightTwoSites(state, references, energy, maxTruncation, deflate, penaltyWeight);

					// right to left
					SweepLeftTwoSites(state, references, energy, maxTruncation, deflate, penaltyWeight);
				}
			}
			else
			{
				for (int sweep = 0; sweep < numSweeps; ++sweep)
				{
					// mixing amplitude: active except on the final sweep (so it converges cleanly)
					const double alpha = (useMixing && sweep < numSweeps - 1) ? perturbationFactor : 0.;

					// left to right
					SweepRightOneSite(state, references, energy, maxTruncation, deflate, penaltyWeight, alpha);

					// right to left
					SweepLeftOneSite(state, references, energy, maxTruncation, deflate, penaltyWeight, alpha);
				}
			}

			truncationError = maxTruncation;

			return energy;
		}

		void MPSDMRGAlgorithm::SweepRightTwoSites(MatrixProductState& state, const std::vector<const MatrixProductState*>& references, double& energy, double& maxTruncation, bool deflate, double penaltyWeight)
		{
			const int chi = static_cast<int>(maxStates);
			const int R = deflate ? static_cast<int>(references.size()) : 0;

			for (int i = 0; i < state.N - 1; ++i)
			{
				const Tensor4 theta = SolvePair(i, state, references, energy, deflate, penaltyWeight);
				maxTruncation = std::max(maxTruncation, state.SplitTwoSite(i, theta, chi, svdThreshold, true));
				Lenv[i + 1] = mpo.UpdateLeft(Lenv[i], state.GetSiteTensor(i));
				for (int r = 0; r < R; ++r)
					OLenv[r][i + 1] = UpdateLeftOverlap(OLenv[r][i], state.GetSiteTensor(i), references[r]->GetSiteTensor(i));
			}
		}

		void MPSDMRGAlgorithm::SweepLeftTwoSites(MatrixProductState& state, const std::vector<const MatrixProductState*>& references, double& energy, double& maxTruncation, bool deflate, double penaltyWeight)
		{
			const int chi = static_cast<int>(maxStates);
			const int R = deflate ? static_cast<int>(references.size()) : 0;

			for (int i = state.N - 2; i >= 0; --i)
			{
				const Tensor4 theta = SolvePair(i, state, references, energy, deflate, penaltyWeight);
				maxTruncation = std::max(maxTruncation, state.SplitTwoSite(i, theta, chi, svdThreshold, false));
				Renv[i] = mpo.UpdateRight(Renv[i + 1], state.GetSiteTensor(i + 1));
				for (int r = 0; r < R; ++r)
					ORenv[r][i] = UpdateRightOverlap(ORenv[r][i + 1], state.GetSiteTensor(i + 1), references[r]->GetSiteTensor(i + 1));
			}
		}

		void MPSDMRGAlgorithm::SweepRightOneSite(MatrixProductState& state, const std::vector<const MatrixProductState*>& references, double& energy, double& maxTruncation, bool deflate, double penaltyWeight, double alpha)
		{
			const int chi = static_cast<int>(maxStates);
			const int R = deflate ? static_cast<int>(references.size()) : 0;

			for (int i = 0; i < state.N; ++i)
			{
				SolveSite(i, state, references, energy, deflate, penaltyWeight);

				if (i < state.N - 1)
				{
					const Eigen::MatrixXd mixing = (alpha > 0.) ? LeftMixing(mpo, Lenv[i], state.GetSiteTensor(i), alpha) : Eigen::MatrixXd();
					maxTruncation = std::max(maxTruncation, state.LeftNormalizeExpanded(i, chi, svdThreshold, mixing));
					Lenv[i + 1] = mpo.UpdateLeft(Lenv[i], state.GetSiteTensor(i));
					for (int r = 0; r < R; ++r)
						OLenv[r][i + 1] = UpdateLeftOverlap(OLenv[r][i], state.GetSiteTensor(i), references[r]->GetSiteTensor(i));
				}
			}
		}

		void MPSDMRGAlgorithm::SweepLeftOneSite(MatrixProductState& state, const std::vector<const MatrixProductState*>& references, double& energy, double& maxTruncation, bool deflate, double penaltyWeight, double alpha)
		{
			const int chi = static_cast<int>(maxStates);
			const int R = deflate ? static_cast<int>(references.size()) : 0;

			for (int i = state.N - 1; i >= 0; --i)
			{
				SolveSite(i, state, references, energy, deflate, penaltyWeight);

				if (i > 0)
				{
					const Eigen::MatrixXd mixing = (alpha > 0.) ? RightMixing(mpo, Renv[i], state.GetSiteTensor(i), alpha) : Eigen::MatrixXd();
					maxTruncation = std::max(maxTruncation, state.RightNormalizeExpanded(i, chi, svdThreshold, mixing));
					Renv[i - 1] = mpo.UpdateRight(Renv[i], state.GetSiteTensor(i));
					for (int r = 0; r < R; ++r)
						ORenv[r][i - 1] = UpdateRightOverlap(ORenv[r][i], state.GetSiteTensor(i), references[r]->GetSiteTensor(i));
				}
			}
		}

		Tensor4 MPSDMRGAlgorithm::SolvePair(int i, MatrixProductState& state, const std::vector<const MatrixProductState*>& references, double& energy, bool deflate, double penaltyWeight)
		{
			const int Dl = state.LeftBondDim(i);
			const int Dr = state.RightBondDim(i + 1);

			TwoSiteEffectiveHamiltonian H;
			H.mpo = &mpo;
			H.L = &Lenv[i];
			H.R = &Renv[i + 1];
			H.d = m_d;
			H.Dl = Dl;
			H.Dr = Dr;

			if (deflate)
			{
				const Eigen::array<Eigen::IndexPair<int>, 1> c1{ Eigen::IndexPair<int>(1, 0) };
				const Eigen::array<Eigen::IndexPair<int>, 1> c2{ Eigen::IndexPair<int>(2, 1) };
				const Eigen::array<Eigen::IndexPair<int>, 1> c3{ Eigen::IndexPair<int>(2, 0) };
				
				std::vector<Eigen::VectorXd> vecs;
				for (int r = 0; r < static_cast<int>(references.size()); ++r)
				{
					const Tensor3& si = references[r]->GetSiteTensor(i);
					const Tensor3& si1 = references[r]->GetSiteTensor(i + 1);

					const Tensor4 u = OLenv[r][i].contract(si, c1).contract(si1.contract(ORenv[r][i + 1], c2), c3);
				    
					assert(u.dimension(0) == Dl && u.dimension(1) == m_d && u.dimension(2) == m_d && u.dimension(3) == Dr);

					Eigen::VectorXd v = Eigen::Map<const Eigen::VectorXd>(u.data(), H.Size());
					const double nrm = v.norm();
					if (nrm > LANCZOS_RESIDUAL_TOLERANCE)
						vecs.push_back(v / nrm);
				}
				if (!vecs.empty())
				{
					H.penaltyVectors.resize(vecs.front().size(), vecs.size());
					for (size_t c = 0; c < vecs.size(); ++c)
						H.penaltyVectors.col(c) = vecs[c];
					H.penaltyWeight = penaltyWeight;
				}
			}

			Eigen::VectorXd guess = H.Pack(state.ContractTwoSite(i));
			energy = LanczosGroundState(H, guess);
			Tensor4 theta = H.Unpack(guess);

			return theta;
		}

		void MPSDMRGAlgorithm::SolveSite(int i, MatrixProductState& state, const std::vector<const MatrixProductState*>& references, double& energy, bool deflate, double penaltyWeight)
		{
			const int Dl = state.LeftBondDim(i);
			const int Dr = state.RightBondDim(i);

			SingleSiteEffectiveHamiltonian H;
			H.mpo = &mpo;
			H.L = &Lenv[i];
			H.R = &Renv[i];
			H.d = m_d;
			H.Dl = Dl;
			H.Dr = Dr;

			if (deflate)
			{
				H.penaltyVectors = BuildSitePenalty(i, Dl, Dr, references, deflate);
				if (H.penaltyVectors.cols() > 0)
					H.penaltyWeight = penaltyWeight;
			}

			Eigen::VectorXd guess = H.Pack(state.GetSiteTensor(i));
			energy = LanczosGroundState(H, guess);
			state.GetSiteTensor(i) = H.Unpack(guess);
		}

		Eigen::MatrixXd MPSDMRGAlgorithm::BuildSitePenalty(int i, int Dl, int Dr, const std::vector<const MatrixProductState*>& references, bool deflate) const
		{
			// Assemble the deflation penalty columns for a single-site tensor at site i.
			Eigen::MatrixXd cols;
			if (!deflate) return cols;

			const int R = deflate ? static_cast<int>(references.size()) : 0;

			const Eigen::array<Eigen::IndexPair<int>, 1> c1{ Eigen::IndexPair<int>(1, 0) };
			const Eigen::array<Eigen::IndexPair<int>, 1> c2{ Eigen::IndexPair<int>(2, 1) };

			std::vector<Eigen::VectorXd> vecs;
			for (int r = 0; r < R; ++r)
			{
				const Tensor3 u = OLenv[r][i].contract(references[r]->GetSiteTensor(i), c1).contract(ORenv[r][i], c2);

				assert (u.dimension(0) == Dl && u.dimension(1) == m_d && u.dimension(2) == Dr);

				Eigen::VectorXd v = Eigen::Map<const Eigen::VectorXd>(u.data(), static_cast<Eigen::Index>(m_d) * Dl * Dr);
				const double nrm = v.norm();
				if (nrm > LANCZOS_RESIDUAL_TOLERANCE)
					vecs.push_back(v / nrm);
			}

			if (!vecs.empty())
			{
				cols.resize(vecs.front().size(), vecs.size());
				for (size_t c = 0; c < vecs.size(); ++c)
					cols.col(c) = vecs[c];
			}

			return cols;
		}

		double MPSDMRGAlgorithm::EnergyExpectation(const MatrixProductState& state)
		{
			const int N = state.N;
			if (N <= 0) return std::numeric_limits<double>::infinity();

			InitREnv(state);

			SingleSiteEffectiveHamiltonian H;
			H.mpo = &mpo;
			Tensor3 L = mpo.LeftBoundary();
			H.L = &L;
			H.R = &Renv[0];
			H.d = m_d;
			H.Dl = state.LeftBondDim(0);
			H.Dr = state.RightBondDim(0);

			const Eigen::VectorXd center = H.Pack(state.GetSiteTensor(0));
			const double norm = center.squaredNorm();
			if (norm <= 0.) return std::numeric_limits<double>::infinity();

			return center.dot(H.Apply(center)) / norm;
		}


		void MPSDMRGAlgorithm::ComputeMeasurements(MatrixProductState& state)
		{
			// NOTE: this sweeps the orthogonality center to the right, so on return
			// 'state' is left-canonical (center at site N-1) rather than right-
			// canonical. That is fine for using it as a deflation reference later,
			// because the deflation only ever reads it through gauge-independent
			// overlaps.
			const int N = state.N;
			const int chi = static_cast<int>(maxStates);

			const Eigen::array<int, 2> transpose{ 1, 0 };
			const Tensor2 Sminus = m_Splus.shuffle(transpose);

			// the state comes in right-canonical form (center at site 0); sweep the
			// orthogonality center to the right, computing the nearest-neighbour bond
			// energy <h_i,i+1> at each bond (this is what the chart plots)
			for (int i = 0; i < N - 1; ++i)
			{
				const double SzSz = state.TwoSiteExpectation(i, m_Sz, m_Sz);
				const double SpSm = state.TwoSiteExpectation(i, m_Splus, Sminus);
				const double SmSp = state.TwoSiteExpectation(i, Sminus, m_Splus);

				results.push_back(m_Jz * SzSz + 0.5 * m_Jxy * (SpSm + SmSp));

				state.LeftNormalize(i, chi, svdThreshold);
			}
		}

	}

}
