#pragma once

#undef min
#undef max
#include <Eigen\eigen>
#include <unsupported/Eigen/CXX11/Tensor>

namespace DMRG {

	namespace Operators {

		using Tensor2 = Eigen::Tensor<double, 2>;

		class Operator
		{
		protected:
			bool changeSign; // not used yet, but will be used for fermionic operators when/if implemented
			int singleSiteSize;
		public:
			// extendChangeSign should be true for fermionic operators (not implemented yet)
			// false for bosonic operators
			Operator(unsigned int size = 2, bool extendChangeSign = false)
				: changeSign(extendChangeSign), singleSiteSize(size)
			{
				matrix = Tensor2(size, size);
				matrix.setZero();
			}

			virtual ~Operator() = default;

			// The operator is stored as a rank-2 tensor so it can be contracted
			// directly with the MPS/MPO tensors without converting to a matrix.
			Tensor2 matrix;

			// Zero-copy Eigen::Map views of the operator storage as a dense matrix.
			// Eigen::Tensor is column-major and contiguous, so the map aliases the
			// tensor data (no copy) and can be used with the usual matrix API
			// (topLeftCorner, block, ...). The tensor must outlive the view.
			Eigen::Map<Eigen::MatrixXd> AsMatrix()
			{
				return Eigen::Map<Eigen::MatrixXd>(matrix.data(), matrix.dimension(0), matrix.dimension(1));
			}

			Eigen::Map<const Eigen::MatrixXd> AsMatrix() const
			{
				return Eigen::Map<const Eigen::MatrixXd>(matrix.data(), matrix.dimension(0), matrix.dimension(1));
			}

			int GetSingleSiteSize() const { return singleSiteSize; }
		};

	}

}