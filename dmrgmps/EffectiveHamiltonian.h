#pragma once

#include "MPO.h"

namespace DMRG {

	namespace MPS {

		class EffectiveHamiltonian
		{
		public:
			const HeisenbergMPO* mpo = nullptr;
			const Tensor3* L = nullptr; // left environment,  L(i, a, j) is (Dl x Dw x Dl)
			const Tensor3* R = nullptr; // right environment, R(i, b, j) is (Dr x Dw x Dr)

			int d = 2;
			int Dl = 1;
			int Dr = 1;

			// deflation: each column of 'penaltyVectors' is a normalized projection of
			// a previously found (reference) state onto the current one-site basis. The
			// term  weight * sum_k |u_k><u_k|  pushes the search away from all of them,
			// which lets several excited states be targeted at once.
			double penaltyWeight = 0.;
			Eigen::MatrixXd penaltyVectors;
		};

		// Effective single site Hamiltonian: acts on a site tensor packed as a flat
		// vector (physical index major, then row-major over the Dl x Dr block).
		// Optionally adds a rank-one deflation penalty  weight * |u><u|  used to
		// target excited states orthogonal to a previously found state.
		class SingleSiteEffectiveHamiltonian : public EffectiveHamiltonian
		{
		public:
			Eigen::Index Size() const { return static_cast<Eigen::Index>(d) * Dl * Dr; }

			//   out(l', sp, r') = sum_{a,b,l,r,s} W(a,b,sp,s) L(l,a,l') R(r,b,r') psi(l,s,r) 
			Eigen::VectorXd Apply(const Eigen::VectorXd& x) const
			{
				// zero-copy view of the flat vector as a (Dl, d, Dr) tensor (both use
				// the same column-major ordering, see Pack/Unpack)
				const Eigen::TensorMap<const Tensor3> psi(x.data(), Dl, d, Dr);

				// contract L's ket index (0) with psi's left index (0)
				const Eigen::array<Eigen::IndexPair<int>, 1> c1{ Eigen::IndexPair<int>(0, 0) };
				const Eigen::Tensor<double, 4> T1 = L->contract(psi, c1); // (a, l', s, r)

				// contract (a) of T1 with W's left MPO index, and psi's physical (input)
				// index s of T1 with W's input index (index 3)
				const Eigen::array<Eigen::IndexPair<int>, 2> c2{ Eigen::IndexPair<int>(0, 0), Eigen::IndexPair<int>(2, 3) };
				const Eigen::Tensor<double, 4> T2 = T1.contract(mpo->W, c2); // (l', r, b, sp)

				// contract (r, b) of T2 with (i=0, b=1) of R(r, b, r')
				const Eigen::array<Eigen::IndexPair<int>, 2> c3{ Eigen::IndexPair<int>(1, 0), Eigen::IndexPair<int>(2, 1) };
				const Tensor3 out = T2.contract(*R, c3); // (l', sp, r')

				// zero-copy view of the tensor storage as the flat output vector
				Eigen::VectorXd result = Eigen::Map<const Eigen::VectorXd>(out.data(), Size());

				if (penaltyWeight != 0. && penaltyVectors.rows() == result.size() && penaltyVectors.cols() > 0)
					result += penaltyWeight * (penaltyVectors * (penaltyVectors.transpose() * x));

				return result;
			}

			// View a flat vector as a rank-3 site tensor psi(l, s, r) without copying.
			// Eigen::Tensor is column-major, so the flat layout is l + Dl*s + Dl*d*r,
			// which is the natural ordering of a (Dl, d, Dr) tensor's storage.
			Eigen::TensorMap<const Tensor3> Unpack(const Eigen::VectorXd& x) const
			{
				return Eigen::TensorMap<const Tensor3>(x.data(), Dl, d, Dr);
			}

			// Copy a rank-3 site tensor psi(l, s, r) into a flat vector (column-major).
			Eigen::VectorXd Pack(const Tensor3& psi) const
			{
				return Eigen::Map<const Eigen::VectorXd>(psi.data(), Size());
			}
		};

		// Effective two-site Hamiltonian used by two-site DMRG. It acts on a two-site
		// wavefunction psi(l, s1, s2, r) of shape (Dl, d, d, Dr), packed column-major.
		// The same bulk MPO tensor W is used on both sites, tied by the shared MPO
		// bond index. Optional multi-vector deflation, as in SingleSiteEffectiveHamiltonian.
		class TwoSiteEffectiveHamiltonian : public EffectiveHamiltonian
		{
		public:
			Eigen::Index Size() const { return static_cast<Eigen::Index>(d) * d * Dl * Dr; }

			//   out(l', s1', s2', r') = sum L(l,a,l') W(a,b,s1',s1) W(b,c,s2',s2)
			//                               psi(l,s1,s2,r) R(r,c,r')
			Eigen::VectorXd Apply(const Eigen::VectorXd& x) const
			{
				const Eigen::TensorMap<const Tensor4> psi(x.data(), Dl, d, d, Dr);

				// contract L's ket index (0) with psi's left index (0)
				const Eigen::array<Eigen::IndexPair<int>, 1> c1{ Eigen::IndexPair<int>(0, 0) };
				const Eigen::Tensor<double, 5> T1 = L->contract(psi, c1); // (a, l', s1, s2, r)

				// contract a (idx0 with W idx0) and s1 (idx2 with W input idx3)
				const Eigen::array<Eigen::IndexPair<int>, 2> c2{ Eigen::IndexPair<int>(0, 0), Eigen::IndexPair<int>(2, 3) };
				const Eigen::Tensor<double, 5> T2 = T1.contract(mpo->W, c2); // (l', s2, r, b, s1')

				// contract b (idx3 with W idx0) and s2 (idx1 with W input idx3)
				const Eigen::array<Eigen::IndexPair<int>, 2> c3{ Eigen::IndexPair<int>(3, 0), Eigen::IndexPair<int>(1, 3) };
				const Eigen::Tensor<double, 5> T3 = T2.contract(mpo->W, c3); // (l', r, s1', c, s2')

				// contract r (idx1 with R idx0) and c (idx3 with R idx1)
				const Eigen::array<Eigen::IndexPair<int>, 2> c4{ Eigen::IndexPair<int>(1, 0), Eigen::IndexPair<int>(3, 1) };
				const Tensor4 out = T3.contract(*R, c4); // (l', s1', s2', r')

				Eigen::VectorXd result = Eigen::Map<const Eigen::VectorXd>(out.data(), Size());

				if (penaltyWeight != 0. && penaltyVectors.rows() == result.size() && penaltyVectors.cols() > 0)
					result += penaltyWeight * (penaltyVectors * (penaltyVectors.transpose() * x));

				return result;
			}

			Eigen::TensorMap<const Tensor4> Unpack(const Eigen::VectorXd& x) const
			{
				return Eigen::TensorMap<const Tensor4>(x.data(), Dl, d, d, Dr);
			}

			Eigen::VectorXd Pack(const Tensor4& psi) const
			{
				return Eigen::Map<const Eigen::VectorXd>(psi.data(), Size());
			}
		};

	}

}
