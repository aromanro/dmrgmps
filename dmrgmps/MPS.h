#pragma once

// Matrix Product State representation used by the MPS-based DMRG engine.
//
// A state on a chain of N sites, each with a 'd'-dimensional physical space, is
// stored as a list of rank-3 site tensors A[site], each an Eigen::Tensor<double,3>
// with dimensions (Dleft, d, Dright) - i.e. A(l, s, r). This mirrors the way the
// QCSim MPS simulator keeps its gamma tensors as Eigen::Tensor<...,3>. The virtual
// (bond) dimension between two neighbour sites is shared.
//
// The class provides the SVD based canonicalization / truncation used while
// sweeping (bond dimension limited by 'chi' and by a singular value threshold,
// the same idea used by the QCSim MPS simulator).

#undef min
#undef max
#include <Eigen/Eigen>
#include <unsupported/Eigen/CXX11/Tensor>

#include <vector>
#include <random>

namespace DMRG {

	namespace MPS {

		using Tensor2 = Eigen::Tensor<double, 2>;
		using Tensor3 = Eigen::Tensor<double, 3>;
		using Tensor4 = Eigen::Tensor<double, 4>;

		class MatrixProductState
		{
		private:
			// chain[site] is a rank-3 tensor with dimensions (Dleft, d, Dright): A(l, s, r)
			std::vector<Tensor3> chain;

		public:
			int N = 0; // number of sites
			int d = 2; // physical dimension

			MatrixProductState() = default;

			MatrixProductState(int nrSites, int physicalDim, int maxBondDim)
			{
				Init(nrSites, physicalDim, maxBondDim);
			}

			MatrixProductState(const MatrixProductState&) = default;
			MatrixProductState(MatrixProductState&&) = default;

			MatrixProductState& operator=(const MatrixProductState&) = default;
			MatrixProductState& operator=(MatrixProductState&&) = default;

			const Tensor3& GetSiteTensor(int site) const
			{
				return chain[site];
			}

			Tensor3& GetSiteTensor(int site)
			{
				return chain[site];
			}

			// Zero-copy Eigen::Map view of the (Dleft x Dright) matrix slice A(:, s, :).
			//
			// Eigen::Tensor is column-major, so the linear offset of A(l, s, r) is
			// l + Dl*s + Dl*d*r. For a fixed physical index 's' the elements are
			// contiguous in 'l' (inner stride 1) and separated by Dl*d between columns
			// (outer stride), which is exactly an Eigen::Map with an OuterStride. No
			// data is copied - the returned map aliases the tensor storage, so the
			// tensor must outlive any use of the returned view.
			using ConstSliceMap = Eigen::Map<const Eigen::MatrixXd, Eigen::Unaligned, Eigen::OuterStride<>>;
			using SliceMap = Eigen::Map<Eigen::MatrixXd, Eigen::Unaligned, Eigen::OuterStride<>>;

			static ConstSliceMap Slice(const Tensor3& T, int s)
			{
				const Eigen::Index Dl = T.dimension(0);
				const Eigen::Index d = T.dimension(1);
				const Eigen::Index Dr = T.dimension(2);
				return ConstSliceMap(T.data() + Dl * s, Dl, Dr, Eigen::OuterStride<>(Dl * d));
			}

			// Zero-copy writable view of the slice A(:, s, :).
			static SliceMap Slice(Tensor3& T, int s)
			{
				const Eigen::Index Dl = T.dimension(0);
				const Eigen::Index d = T.dimension(1);
				const Eigen::Index Dr = T.dimension(2);
				return SliceMap(T.data() + Dl * s, Dl, Dr, Eigen::OuterStride<>(Dl * d));
			}

			// Store a (Dleft x Dright) matrix (or expression) into slice A(:, s, :).
			template <typename Derived>
			static void SetSlice(Tensor3& T, int s, const Eigen::MatrixBase<Derived>& M)
			{
				Slice(T, s) = M;
			}

			void Init(int nrSites, int physicalDim, int maxBondDim, unsigned int seed = 42)
			{
				N = nrSites;
				d = physicalDim;

				chain.assign(N, Tensor3());

				std::mt19937 gen(seed);
				std::uniform_real_distribution<double> dist(-1., 1.);

				int Dleft = 1;
				for (int site = 0; site < N; ++site)
				{
					// bond dimension is limited both by the maximum allowed value and by
					// the exact dimension of the physical space swept so far
					const long long leftExact = ExactBondDim(site, maxBondDim);
					const long long rightExact = ExactBondDim(site + 1, maxBondDim);

					Dleft = static_cast<int>(std::min<long long>(leftExact, maxBondDim));
					if (site == 0) Dleft = 1;

					int Dright = static_cast<int>(std::min<long long>(rightExact, maxBondDim));
					if (site == N - 1) Dright = 1;

					Tensor3 T(Dleft, d, Dright);
					for (int l = 0; l < Dleft; ++l)
						for (int s = 0; s < d; ++s)
							for (int r = 0; r < Dright; ++r)
								T(l, s, r) = dist(gen);

					chain[site] = std::move(T);
				}
			}

			int LeftBondDim(int site) const { return static_cast<int>(chain[site].dimension(0)); }
			int RightBondDim(int site) const { return static_cast<int>(chain[site].dimension(2)); }

			// Left-normalize site 'site' with a truncated SVD and push the remainder
			// (S * V^T) into the next site to the right. Returns the discarded weight.
			double LeftNormalize(int site, int chi, double threshold)
			{
				const int Dl = LeftBondDim(site);
				const int Dr = RightBondDim(site);

				// reshape into a (Dl*d) x Dr matrix, row index = s*Dl + l
				Eigen::MatrixXd M(static_cast<Eigen::Index>(Dl) * d, Dr);
				for (int s = 0; s < d; ++s)
					M.block(static_cast<Eigen::Index>(s) * Dl, 0, Dl, Dr) = Slice(chain[site], s);

				Eigen::BDCSVD<Eigen::MatrixXd, Eigen::ComputeThinU | Eigen::ComputeThinV> svd(M);

				const Eigen::VectorXd& S = svd.singularValues();
				const int keep = ChooseKept(S, chi, threshold);
				const double discarded = DiscardedWeight(S, keep);

				const Eigen::MatrixXd U = svd.matrixU().leftCols(keep);
				const Eigen::MatrixXd V = svd.matrixV().leftCols(keep);
				const Eigen::VectorXd s = S.head(keep);

				// new left-normalized tensor at 'site'
				{
					Tensor3 T(Dl, d, keep);
					for (int ps = 0; ps < d; ++ps)
						SetSlice(T, ps, U.block(static_cast<Eigen::Index>(ps) * Dl, 0, Dl, keep));
					chain[site] = std::move(T);
				}

				// carry = S * V^T  (keep x Dr), fold into next site
				if (site + 1 < N)
				{
					const Eigen::MatrixXd carry = s.asDiagonal() * V.transpose();
					const int Dr2 = RightBondDim(site + 1);
					Tensor3 T(keep, d, Dr2);
					for (int ps = 0; ps < d; ++ps)
						SetSlice(T, ps, carry * Slice(chain[site + 1], ps));
					chain[site + 1] = std::move(T);
				}

				return discarded;
			}

			// Right-normalize site 'site' with a truncated SVD and push the remainder
			// (U * S) into the previous site to the left. Returns the discarded weight.
			double RightNormalize(int site, int chi, double threshold)
			{
				const int Dl = LeftBondDim(site);
				const int Dr = RightBondDim(site);

				// reshape into a Dl x (d*Dr) matrix, col index = s*Dr + r
				Eigen::MatrixXd M(Dl, static_cast<Eigen::Index>(Dr) * d);
				for (int s = 0; s < d; ++s)
					M.block(0, static_cast<Eigen::Index>(s) * Dr, Dl, Dr) = Slice(chain[site], s);

				Eigen::BDCSVD<Eigen::MatrixXd, Eigen::ComputeThinU | Eigen::ComputeThinV> svd(M);

				const Eigen::VectorXd& S = svd.singularValues();
				const int keep = ChooseKept(S, chi, threshold);
				const double discarded = DiscardedWeight(S, keep);

				const Eigen::MatrixXd U = svd.matrixU().leftCols(keep);
				const Eigen::MatrixXd V = svd.matrixV().leftCols(keep); // (d*Dr) x keep
				const Eigen::VectorXd s = S.head(keep);

				// new right-normalized tensor at 'site': A_s(l', r) = V(s*Dr + r, l')
				{
					Tensor3 T(keep, d, Dr);
					for (int ps = 0; ps < d; ++ps)
					{
						Eigen::MatrixXd blk = V.block(static_cast<Eigen::Index>(ps) * Dr, 0, Dr, keep); // Dr x keep
						SetSlice(T, ps, blk.transpose()); // keep x Dr
					}
					chain[site] = std::move(T);
				}

				if (site - 1 >= 0)
				{
					const Eigen::MatrixXd carry = U * s.asDiagonal(); // Dl x keep
					const int Dl2 = LeftBondDim(site - 1);
					Tensor3 T(Dl2, d, keep);
					for (int ps = 0; ps < d; ++ps)
						SetSlice(T, ps, Slice(chain[site - 1], ps) * carry);
					chain[site - 1] = std::move(T);
				}

				return discarded;
			}

			// Left-normalize site 'site' with a truncated SVD, but first widen the
			// (Dl*d) x Dr matrix with an extra 'mixing' block of columns (White's
			// density-matrix / subspace-expansion correction). Only the physical
			// (first Dr) rows of V^T are folded into the next site, so the extra
			// columns just seed new bond directions into U; the represented state is
			// changed by O(||mixing||), which is why the caller shrinks the mixing
			// amplitude to zero on the final sweep.
			double LeftNormalizeExpanded(int site, int chi, double threshold, const Eigen::MatrixXd& mixing)
			{
				if (mixing.size() == 0)
					return LeftNormalize(site, chi, threshold);

				const int Dl = LeftBondDim(site);
				const int Dr = RightBondDim(site);
				const Eigen::Index extra = mixing.cols();

				Eigen::MatrixXd M(static_cast<Eigen::Index>(Dl) * d, static_cast<Eigen::Index>(Dr) + extra);
				for (int s = 0; s < d; ++s)
					M.block(static_cast<Eigen::Index>(s) * Dl, 0, Dl, Dr) = Slice(chain[site], s);
				M.block(0, Dr, static_cast<Eigen::Index>(Dl) * d, extra) = mixing;

				Eigen::BDCSVD<Eigen::MatrixXd, Eigen::ComputeThinU | Eigen::ComputeThinV> svd(M);

				const Eigen::VectorXd& S = svd.singularValues();
				const int keep = ChooseKept(S, chi, threshold);
				const double discarded = DiscardedWeight(S, keep);

				const Eigen::MatrixXd U = svd.matrixU().leftCols(keep);
				const Eigen::MatrixXd V = svd.matrixV().leftCols(keep); // (Dr + extra) x keep
				const Eigen::VectorXd s = S.head(keep);

				{
					Tensor3 T(Dl, d, keep);
					for (int ps = 0; ps < d; ++ps)
						SetSlice(T, ps, U.block(static_cast<Eigen::Index>(ps) * Dl, 0, Dl, keep));
					chain[site] = std::move(T);
				}

				if (site + 1 < N)
				{
					// only the physical part of V (first Dr rows) carries the state
					const Eigen::MatrixXd carry = s.asDiagonal() * V.topRows(Dr).transpose(); // keep x Dr
					const int Dr2 = RightBondDim(site + 1);
					Tensor3 T(keep, d, Dr2);
					for (int ps = 0; ps < d; ++ps)
						SetSlice(T, ps, carry * Slice(chain[site + 1], ps));
					chain[site + 1] = std::move(T);
				}

				return discarded;
			}

			// Right-normalizing counterpart of LeftNormalizeExpanded: widen the
			// Dl x (d*Dr) matrix with an extra 'mixing' block of rows and fold only
			// the physical (first Dl) columns of U into the previous site.
			double RightNormalizeExpanded(int site, int chi, double threshold, const Eigen::MatrixXd& mixing)
			{
				if (mixing.size() == 0)
					return RightNormalize(site, chi, threshold);

				const int Dl = LeftBondDim(site);
				const int Dr = RightBondDim(site);
				const Eigen::Index extra = mixing.rows();

				Eigen::MatrixXd M(static_cast<Eigen::Index>(Dl) + extra, static_cast<Eigen::Index>(Dr) * d);
				for (int s = 0; s < d; ++s)
					M.block(0, static_cast<Eigen::Index>(s) * Dr, Dl, Dr) = Slice(chain[site], s);
				M.block(Dl, 0, extra, static_cast<Eigen::Index>(Dr) * d) = mixing;

				Eigen::BDCSVD<Eigen::MatrixXd, Eigen::ComputeThinU | Eigen::ComputeThinV> svd(M);

				const Eigen::VectorXd& S = svd.singularValues();
				const int keep = ChooseKept(S, chi, threshold);
				const double discarded = DiscardedWeight(S, keep);

				const Eigen::MatrixXd U = svd.matrixU().leftCols(keep); // (Dl + extra) x keep
				const Eigen::MatrixXd V = svd.matrixV().leftCols(keep); // (d*Dr) x keep
				const Eigen::VectorXd s = S.head(keep);

				{
					Tensor3 T(keep, d, Dr);
					for (int ps = 0; ps < d; ++ps)
					{
						Eigen::MatrixXd blk = V.block(static_cast<Eigen::Index>(ps) * Dr, 0, Dr, keep);
						SetSlice(T, ps, blk.transpose());
					}
					chain[site] = std::move(T);
				}

				if (site - 1 >= 0)
				{
					// only the physical part of U (first Dl rows) carries the state
					const Eigen::MatrixXd carry = U.topRows(Dl) * s.asDiagonal(); // Dl x keep
					const int Dl2 = LeftBondDim(site - 1);
					Tensor3 T(Dl2, d, keep);
					for (int ps = 0; ps < d; ++ps)
						SetSlice(T, ps, Slice(chain[site - 1], ps) * carry);
					chain[site - 1] = std::move(T);
				}

				return discarded;
			}

			// Contract the neighbour site tensors chain[site] (Dl,d,m) and chain[site+1]
			// (m,d,Dr) over the shared bond into a two-site wavefunction
			// theta(l, s1, s2, r) of shape (Dl, d, d, Dr).
			Tensor4 ContractTwoSite(int site) const
			{
				const Eigen::array<Eigen::IndexPair<int>, 1> c{ Eigen::IndexPair<int>(2, 0) };
				return chain[site].contract(chain[site + 1], c);
			}

			// Split an optimized two-site tensor theta(l, s1, s2, r) back into
			// chain[site] and chain[site+1] with a truncated SVD (bond limited by chi). When
			// 'centerRight' the singular values go to chain[site+1] (chain[site] becomes
			// left-normalized); otherwise they go to chain[site] (chain[site+1] becomes
			// right-normalized). Returns the discarded weight.
			double SplitTwoSite(int site, const Tensor4& theta, int chi, double threshold, bool centerRight)
			{
				const int Dl = static_cast<int>(theta.dimension(0));
				const int Dr = static_cast<int>(theta.dimension(3));

				// reshape into (Dl*d) x (d*Dr): row = s1*Dl + l, col = s2*Dr + r
				Eigen::MatrixXd M(static_cast<Eigen::Index>(Dl) * d, static_cast<Eigen::Index>(Dr) * d);
				for (int s1 = 0; s1 < d; ++s1)
					for (int s2 = 0; s2 < d; ++s2)
						for (int l = 0; l < Dl; ++l)
							for (int r = 0; r < Dr; ++r)
								M(static_cast<Eigen::Index>(s1) * Dl + l, static_cast<Eigen::Index>(s2) * Dr + r) = theta(l, s1, s2, r);

				Eigen::BDCSVD<Eigen::MatrixXd, Eigen::ComputeThinU | Eigen::ComputeThinV> svd(M);

				const Eigen::VectorXd& S = svd.singularValues();
				const int keep = ChooseKept(S, chi, threshold);
				const double discarded = DiscardedWeight(S, keep);

				const Eigen::MatrixXd U = svd.matrixU().leftCols(keep); // (Dl*d) x keep
				const Eigen::MatrixXd V = svd.matrixV().leftCols(keep); // (d*Dr) x keep
				const Eigen::VectorXd s = S.head(keep);

				Tensor3 TL(Dl, d, keep);
				for (int s1 = 0; s1 < d; ++s1)
				{
					const Eigen::MatrixXd blk = U.block(static_cast<Eigen::Index>(s1) * Dl, 0, Dl, keep); // Dl x keep
					SetSlice(TL, s1, centerRight ? blk : (blk * s.asDiagonal()).eval());
				}

				Tensor3 TR(keep, d, Dr);
				for (int s2 = 0; s2 < d; ++s2)
				{
					const Eigen::MatrixXd blk = V.block(static_cast<Eigen::Index>(s2) * Dr, 0, Dr, keep); // Dr x keep
					SetSlice(TR, s2, centerRight ? (s.asDiagonal() * blk.transpose()).eval() : blk.transpose().eval());
				}

				chain[site] = std::move(TL);
				chain[site + 1] = std::move(TR);

				return discarded;
			}

			// Bring the whole chain to right canonical form (orthogonality center at site 0).
			void RightCanonicalize(int chi, double threshold)
			{
				for (int site = N - 1; site > 0; --site)
					RightNormalize(site, chi, threshold);
			}

			// Expectation value <op> for the single-site operator 'op' acting on 'site',
			// assuming the orthogonality center is at 'site' (state normalized).
			double SingleSiteExpectation(int site, const Tensor2& op) const
			{
				std::vector<Eigen::MatrixXd> a(d);
				for (int s = 0; s < d; ++s) a[s] = Slice(chain[site], s);

				double result = 0.;
				double norm = 0.;
				for (int sp = 0; sp < d; ++sp)
				{
					norm += (a[sp].array() * a[sp].array()).sum();
					for (int s = 0; s < d; ++s)
					{
						const double coef = op(sp, s);
						if (coef == 0.) continue;
						result += coef * (a[sp].array() * a[s].array()).sum();
					}
				}

				if (norm <= 0.) return 0.;

				return result / norm;
			}

			// Expectation value <op1_i op2_{i+1}> for operators acting on neighbour
			// sites 'site' and 'site+1', assuming the orthogonality center is at 'site'
			// (so A[site] is the center and A[site+1] is right-normalized).
			double TwoSiteExpectation(int site, const Tensor2& op1, const Tensor2& op2) const
			{
				std::vector<Eigen::MatrixXd> a(d), b(d);
				for (int s = 0; s < d; ++s) { a[s] = Slice(chain[site], s); b[s] = Slice(chain[site + 1], s); }

				// B[sp][sq] = chain[site][sp] * chain[site+1][sq]  (Dl x Dr matrix)
				std::vector<std::vector<Eigen::MatrixXd>> B(d, std::vector<Eigen::MatrixXd>(d));
				for (int sp = 0; sp < d; ++sp)
					for (int sq = 0; sq < d; ++sq)
						B[sp][sq] = a[sp] * b[sq];

				double result = 0.;
				double norm = 0.;
				for (int sp = 0; sp < d; ++sp)
					for (int sq = 0; sq < d; ++sq)
						norm += (B[sp][sq].array() * B[sp][sq].array()).sum();

				for (int sp = 0; sp < d; ++sp)
					for (int sq = 0; sq < d; ++sq)
					{
						const double c1 = op1(sp, sq);
						for (int tp = 0; tp < d; ++tp)
							for (int tq = 0; tq < d; ++tq)
							{
								const double c2 = op2(tp, tq);
								const double coef = c1 * c2;
								if (coef == 0.) continue;
								// <B[sp][tp], B[sq][tq]> Frobenius inner product
								result += coef * (B[sp][tp].array() * B[sq][tq].array()).sum();
							}
					}

				if (norm <= 0.) return 0.;

				return result / norm;
			}

		private:
			long long ExactBondDim(int cut, long long maxBondDim) const
			{
				// dimension of the smaller of the two halves at bond 'cut'
				long long left = 1;
				for (int i = 0; i < cut && left <= maxBondDim; ++i) left *= d;
				long long right = 1;
				for (int i = cut; i < N && right <= maxBondDim; ++i) right *= d;

				return std::min<long long>(left, right);
			}

			static int ChooseKept(const Eigen::VectorXd& S, int chi, double threshold)
			{
				int nonzero = 0;
				const double maxS = (S.size() > 0) ? S(0) : 0.;
				const double cut = threshold * maxS;
				for (Eigen::Index i = 0; i < S.size(); ++i)
					if (S(i) > cut && S(i) > 0.) ++nonzero;

				if (nonzero < 1) nonzero = 1;

				int keep = nonzero;
				if (chi > 0 && keep > chi) keep = chi;

				return keep;
			}

			static double DiscardedWeight(const Eigen::VectorXd& S, int keep)
			{
				double total = 0.;
				double kept = 0.;
				for (Eigen::Index i = 0; i < S.size(); ++i)
				{
					const double w = S(i) * S(i);
					total += w;
					if (i < keep) kept += w;
				}

				if (total <= 0.) return 0.;

				return 1. - kept / total;
			}
		};

	}

}
