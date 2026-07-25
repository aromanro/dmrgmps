#pragma once

#include "Operator.h"

#include <stdexcept>


namespace DMRG {

	namespace Operators {

		inline unsigned int ValidateSpinOneHalfSize(unsigned int size)
		{
			if (size == 0 || size % 2 != 0)
				throw std::invalid_argument("Spin-1/2 operator size must be a positive multiple of 2.");

			return size;
		}

		inline unsigned int ValidateSpinOneSize(unsigned int size)
		{
			if (size == 0 || size % 3 != 0)
				throw std::invalid_argument("Spin-1 operator size must be a positive multiple of 3.");

			return size;
		}

		class SiteOperator :
			public Operator
		{
		public:
			SiteOperator(unsigned int size = 2, bool extendChangeSign = false)
				: Operator(size, extendChangeSign)
			{}
		};


		class SzOneHalf : public SiteOperator
		{
		public:
			SzOneHalf(unsigned int size = 2)
				: SiteOperator(ValidateSpinOneHalfSize(size), false)
			{
				const int subsize = size / 2;

				Eigen::Map<Eigen::MatrixXd> m = AsMatrix();
				m.topLeftCorner(subsize, subsize) = 0.5 * Eigen::MatrixXd::Identity(subsize, subsize);
				m.bottomRightCorner(subsize, subsize) = -0.5 * Eigen::MatrixXd::Identity(subsize, subsize);
			}
		};


		class SplusOneHalf : public SiteOperator
		{
		public:
			SplusOneHalf(unsigned int size = 2)
				: SiteOperator(ValidateSpinOneHalfSize(size), false)
			{
				const int subsize = size / 2;

				AsMatrix().topRightCorner(subsize, subsize) = Eigen::MatrixXd::Identity(subsize, subsize);
			}
		};


		class SzOne : public SiteOperator
		{
		public:
			SzOne(unsigned int size = 3)
				: SiteOperator(ValidateSpinOneSize(size), false)
			{
				const int subsize = size / 3;

				Eigen::Map<Eigen::MatrixXd> m = AsMatrix();
				m.topLeftCorner(subsize, subsize) = Eigen::MatrixXd::Identity(subsize, subsize);
				m.bottomRightCorner(subsize, subsize) = -Eigen::MatrixXd::Identity(subsize, subsize);
			}
		};


		class SplusOne : public SiteOperator
		{
		public:
			SplusOne(unsigned int size = 3)
				: SiteOperator(ValidateSpinOneSize(size), false)
			{
				const int subsize = size / 3;

				// In the {|+1>, |0>, |-1>} basis S+ acts on the superdiagonal blocks:
				// S+|0> = sqrt(2)|+1>, S+|-1> = sqrt(2)|0>.
				Eigen::Map<Eigen::MatrixXd> m = AsMatrix();
				m.block(0, subsize, subsize, subsize) = Eigen::MatrixXd::Identity(subsize, subsize);
				m.block(subsize, 2 * subsize, subsize, subsize) = Eigen::MatrixXd::Identity(subsize, subsize);

				m *= sqrt(2.);
			}
		};

	}
}
