#pragma once

// Matrix Product Operator (MPO) for the nearest-neighbour XXZ Heisenberg chain,
// together with the left/right environment tensors and the effective, single-site
// Hamiltonian used by the variational MPS-DMRG sweeps.
//
// This uses Eigen::Tensor and Eigen::Tensor::contract for the contractions, in the
// same spirit as the QCSim MPS simulator (which keeps its site tensors as
// Eigen::Tensor<...,3> and contracts them with Eigen::Tensor::contract).
//
// Conventions:
//   - MPS site tensor A has dimensions (Dleft, d, Dright): A(l, s, r).
//   - The bulk MPO is a rank-4 tensor W(a, b, sp, s), where a/b are the left/right
//     MPO bond indices (dimension Dw = 5) and sp/s are the output/input physical
//     indices.
//   - Environments are rank-3 tensors: the left environment L(i, a, j) and the
//     right environment R(i, b, j), where i is the ket virtual index, j the bra
//     virtual index and a/b the MPO bond index.
//
// The Hamiltonian reproduced here is exactly the one used by the old block DMRG:
//   H = sum_i  Jz * Sz_i Sz_{i+1} + (Jxy/2) * ( S+_i S-_{i+1} + S-_i S+_{i+1} )
//
// The bulk MPO tensor has bond dimension 5:
//   W(0,0)=I
//   W(1,0)=S+      
//   W(2,0)=S-      
//   W(3,0)=Sz
//   W(4,0)=-h*Sz   W(4,1)=(Jxy/2)S-   W(4,2)=(Jxy/2)S+   W(4,3)=Jz*Sz   W(4,4)=I
//	 See: The density-matrix renormalization group in the age of matrix product states
//   by Ulrich Schollwock
//   6.1. MPO representation of Hamiltonians, formula 184

#undef min
#undef max
#include <Eigen/Eigen>
#include <unsupported/Eigen/CXX11/Tensor>

#include <vector>

#include "MPS.h"

namespace DMRG {

	namespace MPS {

		class HeisenbergMPO
		{
		public:
			static constexpr int Dw = 5;

			int d = 2;

			// W(a, b, sp, s): rank-4 bulk MPO tensor
			// dimensions: (Dw, Dw, d, d)
			Tensor4 W;

			HeisenbergMPO() = default;

			void Build(int physicalDim, const Tensor2& Sz, const Tensor2& Splus, double Jz, double Jxy, double h = 0.)
			{
				d = physicalDim;

				const Eigen::array<int, 2> transpose{ 1, 0 };
				const Tensor2 Sminus = Splus.shuffle(transpose);
				Tensor2 Id(d, d);
				Id.setZero();
				for (int i = 0; i < d; ++i) Id(i, i) = 1.;

				W = Tensor4(Dw, Dw, d, d);
				W.setZero();

				SetBlock(0, 0, Id);
				SetBlock(1, 0, Splus);
				SetBlock(2, 0, Sminus);
				SetBlock(3, 0, Sz);

				SetBlock(4, 0, Tensor2(-h * Sz));

				SetBlock(4, 1, Tensor2(0.5 * Jxy * Sminus));
				SetBlock(4, 2, Tensor2(0.5 * Jxy * Splus));
				SetBlock(4, 3, Tensor2(Jz * Sz));
				SetBlock(4, 4, Id);
			}

			// Left boundary environment (before site 0): bond dimensions 1, MPO index
			// selects the last MPO row.
			Tensor3 LeftBoundary() const
			{
				Tensor3 L(1, Dw, 1);
				L.setZero();
				L(0, Dw - 1, 0) = 1.;
				return L;
			}

			// Right boundary environment (after the last site): bond dimensions 1, MPO
			// index selects the first MPO column.
			Tensor3 RightBoundary() const
			{
				Tensor3 R(1, Dw, 1);
				R.setZero();
				R(0, 0, 0) = 1.;
				return R;
			}

			// Grow the left environment by absorbing a (left-normalized) site tensor A.
			// L(i, a, j) is (Dl x Dw x Dl); result Lnew(i', b, j') is (Dr x Dw x Dr).
			//   Lnew(i', b, j') = sum_{i,a,j,sp,s} L(i,a,j) W(a,b,sp,s) A(i,sp,i') A(j,s,j')
			Tensor3 UpdateLeft(const Tensor3& L, const Tensor3& A) const
			{
				// contract L's ket index (0) with the ket A's left index (0)
				const Eigen::array<Eigen::IndexPair<int>, 1> c1{ Eigen::IndexPair<int>(0, 0) };
				const Eigen::Tensor<double, 4> T1 = L.contract(A, c1); // (a, j, s, i')

				// contract a (T1 idx0 with W idx0) and the ket physical s (T1 idx2 with
				// W's input index 3)
				const Eigen::array<Eigen::IndexPair<int>, 2> c2{ Eigen::IndexPair<int>(0, 0), Eigen::IndexPair<int>(2, 3) };
				const Eigen::Tensor<double, 4> T2 = T1.contract(W, c2); // (j, i', b, sp)

				// contract the bra virtual j (T2 idx0 with A idx0) and the bra physical
				// sp (T2 idx3 with A idx1)
				const Eigen::array<Eigen::IndexPair<int>, 2> c3{ Eigen::IndexPair<int>(0, 0), Eigen::IndexPair<int>(3, 1) };
				Tensor3 Lnew = T2.contract(A, c3); // (i', b, j')

				return Lnew;
			}

			// Grow the right environment by absorbing a (right-normalized) site tensor A.
			// R(i, b, j) is (Dr x Dw x Dr); result Rnew(i', a, j') is (Dl x Dw x Dl).
			//   Rnew(i', a, j') = sum_{i,b,j,sp,s} R(i,b,j) W(a,b,sp,s) A(i',sp,i) A(j',s,j)
			Tensor3 UpdateRight(const Tensor3& R, const Tensor3& A) const
			{
				// contract R's ket index (0) with the ket A's right index (2)
				const Eigen::array<Eigen::IndexPair<int>, 1> c1{ Eigen::IndexPair<int>(0, 2) };
				const Eigen::Tensor<double, 4> T1 = R.contract(A, c1); // (b, j, i', s)

				// contract b (T1 idx0 with W idx1) and the ket physical s (T1 idx3 with
				// W's input index 3)
				const Eigen::array<Eigen::IndexPair<int>, 2> c2{ Eigen::IndexPair<int>(0, 1), Eigen::IndexPair<int>(3, 3) };
				const Eigen::Tensor<double, 4> T2 = T1.contract(W, c2); // (j, i', a, sp)

				// contract the bra virtual j (T2 idx0 with A idx2) and the bra physical
				// sp (T2 idx3 with A idx1)
				const Eigen::array<Eigen::IndexPair<int>, 2> c3{ Eigen::IndexPair<int>(0, 2), Eigen::IndexPair<int>(3, 1) };
				Tensor3 Rnew = T2.contract(A, c3); // (i', a, j')

				return Rnew;
			}

		private:
			void SetBlock(int a, int b, const Tensor2& op)
			{
				for (int sp = 0; sp < d; ++sp)
					for (int s = 0; s < d; ++s)
						W(a, b, sp, s) = op(sp, s);
			}
		};

	}

}

