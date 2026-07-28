// Console sanity check for the MPS-DMRG engine.
//
// This is a non-GUI harness that links directly against the MPS-DMRG core
// (MPSDMRGAlgorithm.cpp / MPS.h / MPO.h). It runs a couple of well known
// Heisenberg benchmarks and checks the ground state energy per bond against
// the accepted reference values:
//
//   * spin-1/2 antiferromagnetic Heisenberg chain (Bethe ansatz, thermodynamic
//     limit): e0 = 1/4 - ln(2) ~= -0.4431471805599...  (energy per bond)
//   * spin-1 antiferromagnetic Heisenberg chain (Haldane): e0 ~= -1.401484039
//     (energy per bond)
//
// For a finite open chain the per-bond energy differs from the thermodynamic
// limit, so the check uses tolerant bounds and mainly verifies that the solver
// converges to a sensible variational energy that decreases with bond dimension.

#include "MPSDMRGAlgorithm.h"

#include <Eigen/Eigen>

#define _USE_MATH_DEFINES 1
#include <math.h>
#include <iomanip>
#include <iostream>
#include <string>

namespace {

	// Spin-1/2 operators in the {|up>, |down>} basis.
	Eigen::MatrixXd SzOneHalf()
	{
		Eigen::MatrixXd m = Eigen::MatrixXd::Zero(2, 2);
		m(0, 0) = 0.5;
		m(1, 1) = -0.5;
		return m;
	}

	Eigen::MatrixXd SplusOneHalf()
	{
		Eigen::MatrixXd m = Eigen::MatrixXd::Zero(2, 2);
		m(0, 1) = 1.;
		return m;
	}

	// Spin-1 operators in the {|+1>, |0>, |-1>} basis, matching SiteOperator.cpp.
	Eigen::MatrixXd SzOne()
	{
		Eigen::MatrixXd m = Eigen::MatrixXd::Zero(3, 3);
		m(0, 0) = 1.;
		m(2, 2) = -1.;
		return m;
	}

	Eigen::MatrixXd SplusOne()
	{
		// {|+1>, |0>, |-1>} basis: S+ lives on the superdiagonal.
		Eigen::MatrixXd m = Eigen::MatrixXd::Zero(3, 3);
		m(0, 1) = 1.;
		m(1, 2) = 1.;
		m *= std::sqrt(2.);
		return m;
	}

	struct CheckResult
	{
		bool passed = true;
		int failures = 0;
	};

	void Check(CheckResult& cr, bool condition, const std::string& message)
	{
		std::cout << (condition ? "  [PASS] " : "  [FAIL] ") << message << '\n';
		if (!condition)
		{
			cr.passed = false;
			++cr.failures;
		}
	}

	// Runs a finite-chain calculation and returns the total ground state energy.
	// results now hold the per-bond correlation <S_i . S_{i+1}> (what the chart
	// plots). For the isotropic (Jz = Jxy = 1) Heisenberg model the sum of the
	// bond correlations equals the total energy, which is a strong consistency
	// check between the reported energy and the plotted data.
	// Convert a dense matrix into a rank-2 tensor (both are column-major, so this
	// is a straight copy of the storage into the tensor the engine now expects).
	Eigen::Tensor<double, 2> ToTensor(const Eigen::MatrixXd& m)
	{
		Eigen::Tensor<double, 2> t(m.rows(), m.cols());
		Eigen::Map<Eigen::MatrixXd>(t.data(), m.rows(), m.cols()) = m;
		return t;
	}

	double RunChain(int physicalDim, const Eigen::MatrixXd& Sz, const Eigen::MatrixXd& Splus,
		int sites, int sweeps, unsigned int maxStates, double Jz, double Jxy,
		double& energyPerBond, double& bondCorrSum, int& bondCount, int method = 0, std::vector<std::vector<double>> *correlations = nullptr)
	{
		DMRG::MPS::MPSDMRGAlgorithm dmrg(physicalDim, ToTensor(Sz), ToTensor(Splus), Jz, Jxy, maxStates, 0, method);

		if (correlations)
		{
			// Add a correlation measurement for <S+_i S-_j> between all bonds.
			const Eigen::MatrixXd Sminus = Splus.transpose();
			dmrg.correlations.emplace_back(ToTensor(Splus), ToTensor(Sminus));
		}

		const double energy = dmrg.CalculateFinite(sites, sweeps);

		energyPerBond = energy / (sites - 1);

		bondCorrSum = 0.;
		for (double s : dmrg.results)
			bondCorrSum += s;

		bondCount = static_cast<int>(dmrg.results.size());

		if (correlations)
			correlations->swap(dmrg.correlationResults);

		return energy;
	}

} // anonymous namespace

int main()
{
	std::cout << std::fixed << std::setprecision(8);
	std::cout << "==================================================\n";
	std::cout << " MPS-DMRG sanity check / console test\n";
	std::cout << "==================================================\n\n";

	CheckResult cr;

	const Eigen::MatrixXd Sz12 = SzOneHalf();
	const Eigen::MatrixXd Sp12 = SplusOneHalf();
	const Eigen::MatrixXd Sz1 = SzOne();
	const Eigen::MatrixXd Sp1 = SplusOne();

	// ---- Sanity checks on the operator algebra --------------------------------
	std::cout << "-- Operator algebra sanity checks --\n";
	{
		// [Sz, S+] = S+  (S+ raises Sz by one unit)
		const Eigen::MatrixXd comm12 = Sz12 * Sp12 - Sp12 * Sz12;
		Check(cr, (comm12 - Sp12).norm() < 1e-12, "spin-1/2: [Sz, S+] = S+");

		const Eigen::MatrixXd comm1 = Sz1 * Sp1 - Sp1 * Sz1;
		Check(cr, (comm1 - Sp1).norm() < 1e-12, "spin-1: [Sz, S+] = S+");

		// S^2 = Sz^2 + 1/2 (S+ S- + S- S+) should be s(s+1) I.
		const Eigen::MatrixXd Sm12 = Sp12.transpose();
		const Eigen::MatrixXd S2_12 = Sz12 * Sz12 + 0.5 * (Sp12 * Sm12 + Sm12 * Sp12);
		Check(cr, (S2_12 - 0.75 * Eigen::MatrixXd::Identity(2, 2)).norm() < 1e-12,
			"spin-1/2: S^2 = 3/4 I");

		const Eigen::MatrixXd Sm1 = Sp1.transpose();
		const Eigen::MatrixXd S2_1 = Sz1 * Sz1 + 0.5 * (Sp1 * Sm1 + Sm1 * Sp1);
		Check(cr, (S2_1 - 2.0 * Eigen::MatrixXd::Identity(3, 3)).norm() < 1e-12,
			"spin-1: S^2 = 2 I");
	}
	std::cout << '\n';

	// ---- Two-site spin-1/2 chain (exact) --------------------------------------
	// Two spins with H = J (Sx Sx + Sy Sy + Sz Sz): singlet ground energy = -3/4 J.
	std::cout << "-- Two-site spin-1/2 Heisenberg (exact = -0.75) --\n";
	{
		for (int method = 0; method < 3; ++method)
		{
			std::cout << "  Method " << method << ": "
				<< (method == 0 ? "single-site" : (method == 1 ? "one-site with subspace expansion" : "two-site")) << '\n';
			double perBond = 0., bondCorrSum = 0.; int bondCount = 0;
			const double e = RunChain(2, Sz12, Sp12, 2, 8, 16, 1., 1., perBond, bondCorrSum, bondCount, method);
			std::cout << "  E(total) = " << e << ",  sum<S.S> = " << bondCorrSum
				<< ",  bonds = " << bondCount << '\n';
			Check(cr, std::abs(e - (-0.75)) < 1e-5, "two-site energy ~= -0.75");
			Check(cr, bondCount == 1, "two-site chart has 1 bond point");
			Check(cr, std::abs(bondCorrSum - e) < 1e-6, "sum<S.S> == E (chart data matches energy)");
		}
	}
	std::cout << '\n';

	// ---- Longer spin-1/2 chain (Bethe ansatz per-bond limit ~ -0.4431) --------
	std::cout << "-- Spin-1/2 Heisenberg chain convergence --\n";
	{
		const double betheLimit = 0.25 - std::log(2.0); // ~= -0.44314718
		std::cout << "  Bethe ansatz per-bond limit = " << betheLimit << '\n';
		const int sites = 120, sweeps = 8;
		const unsigned int bonds[] = { 8, 16, 32 };

		for (int method = 0; method < 3; ++method)
		{
			std::cout << "  Method " << method << ": "
				<< (method == 0 ? "single-site" : (method == 1 ? "one-site with subspace expansion" : "two-site")) << '\n';

			double prevPerBond = 0.;
			bool haveDecrease = true;

			for (int i = 0; i < 3; ++i)
			{
				double perBond = 0., bondCorrSum = 0.; int bondCount = 0;
				const double e = RunChain(2, Sz12, Sp12, sites, sweeps, bonds[i], 1., 1., perBond, bondCorrSum, bondCount, method);
				std::cout << "  m = " << std::setw(3) << bonds[i]
					<< "  E/bond = " << perBond
					<< "  E(total) = " << e
					<< "  sum<S.S> = " << bondCorrSum << '\n';
				if (i > 0 && perBond > prevPerBond + 1e-6)
					haveDecrease = false;
				prevPerBond = perBond;

				if (i == 2) {
					// A finite open chain sits close to the thermodynamic per-bond limit,
					// but is not bounded by it, so only require it to be in the ballpark.
					Check(cr, std::abs(perBond - betheLimit) < 0.01, "E/bond within 0.01 of Bethe limit");
					Check(cr, bondCount == sites - 1, "chart has N-1 bond points");
					Check(cr, std::abs(bondCorrSum - e) < 1e-4, "sum<S.S> == E (chart data matches energy)");
				}
			}
			Check(cr, haveDecrease, "E/bond is non-increasing with bond dimension");
		}
	}
	std::cout << '\n';

	// ---- Spin-1/2 XY even number of sites chain, with known analytical results --------
	// TODO: Also check s^+_i s^-_{i+1} correlations against the analytical result.

	std::cout << "-- Spin-1/2 XY chain with even number of sites convergence --\n";
	{
		const int sites = 20, sweeps = 8;
		const unsigned int bonds[] = { 8, 16, 32 };

		const double energyAnalyticalValue = 0.5 * (1. - 1. / sin(M_PI_2 / (sites + 1.))); 
		std::cout << "  Analytical energy value for 20 sites = " << energyAnalyticalValue << ", per-bond = " << energyAnalyticalValue / (sites - 1) << '\n';

		// for s+s- it's something like:
		// corrVal = 0.5 * 1. / (sites + 1.) * (1. / sin(M_PI_2 / (sites + 1.)) - pow(-1., l + 1.) / sin(M_PI_2 * (2. * l + 3.) / (sites + 1.)));
		// where l = bond index (0..N-2) for the N-site chain.
		// I need to do some slight changes to RunChain to return the s+s- correlations for each bond
		// then I'll add the check

		for (int method = 0; method < 3; ++method)
		{
			std::cout << "  Method " << method << ": "
				<< (method == 0 ? "single-site" : (method == 1 ? "one-site with subspace expansion" : "two-site")) << '\n';

			double prevPerBond = 0.;
			bool haveDecrease = true;

			for (int i = 0; i < 3; ++i)
			{
				double perBond = 0., bondCorrSum = 0.; int bondCount = 0;
				std::vector<std::vector<double>> correlations;
				const double e = RunChain(2, Sz12, Sp12, sites, sweeps, bonds[i], 0., 1., perBond, bondCorrSum, bondCount, method, &correlations);
				std::cout << "  m = " << std::setw(3) << bonds[i]
					<< "  E/bond = " << perBond
					<< "  E(total) = " << e
					<< "  sum<S.S> = " << bondCorrSum << '\n';
				if (i > 0 && perBond > prevPerBond + 1e-6)
					haveDecrease = false;
				prevPerBond = perBond;

				if (i == 2) {
					// A finite open chain sits close to the thermodynamic per-bond limit,
					// but is not bounded by it, so only require it to be in the ballpark.
					Check(cr, std::abs(e - energyAnalyticalValue) < 0.01, "E within 0.01 of analytical value");
					Check(cr, bondCount == sites - 1, "chart has N-1 bond points");
					Check(cr, std::abs(bondCorrSum - e) < 1e-4, "sum<S.S> == E (chart data matches energy)");

					for (int l = 0; l < bondCount; ++l)
					{
						const double corrVal = -0.5 * 1. / (sites + 1.) * (1. / sin(M_PI_2 / (sites + 1.)) - pow(-1., l + 1.) / sin(M_PI_2 * (2. * l + 3.) / (sites + 1.)));
						std::cout << "    bond " << l << ": <S+_i S-_j> = " << correlations[l][0] << ", analytical = " << corrVal << '\n';
						Check(cr, std::abs(correlations[l][0] - corrVal) < 0.005, "bond correlation matches analytical value");
					}
				}
			}
			Check(cr, haveDecrease, "E/bond is non-increasing with bond dimension");
		}
	}
	std::cout << '\n';

	// ---- Spin-1 Haldane chain -------------------------------------------------
	std::cout << "-- Spin-1 Heisenberg chain (Haldane, per-bond ~ -1.4015) --\n";
	{
		const double haldaneLimit = -1.401484039;
		std::cout << "  Reference per-bond limit = " << haldaneLimit << '\n';
		const int sites = 120, sweeps = 8;

		for (int method = 0; method < 3; ++method)
		{
			std::cout << "  Method " << method << ": "
				<< (method == 0 ? "single-site" : (method == 1 ? "one-site with subspace expansion" : "two-site")) << '\n';
			double perBond = 0., bondCorrSum = 0.; int bondCount = 0;
			const double e = RunChain(3, Sz1, Sp1, sites, sweeps, 32, 1., 1., perBond, bondCorrSum, bondCount, method);
			std::cout << "  m =  32  E/bond = " << perBond
				<< "  E(total) = " << e
				<< "  sum<S.S> = " << bondCorrSum << '\n';

			Check(cr, std::abs(perBond - haldaneLimit) < 0.01, "E/bond within 0.01 of Haldane limit");
			Check(cr, bondCount == sites - 1, "chart has N-1 bond points");
			Check(cr, std::abs(bondCorrSum - e) < 1e-3, "sum<S.S> == E (chart data matches energy)");
		}
	}
	std::cout << '\n';

	// ---- Two-site DMRG vs single-site -----------------------------------------
	// Two-site DMRG should reach an energy at least as low as the single-site
	// update on the same chain (it explores a larger variational space per step).
	std::cout << "-- Two-site DMRG vs single-site (spin-1/2, 20 sites) --\n";
	{
		const int sites = 20, sweeps = 5;
		const unsigned int m = 20;

		DMRG::MPS::MPSDMRGAlgorithm single(2, ToTensor(Sz12), ToTensor(Sp12), 1., 1., m);
		single.sweepMode = DMRG::MPS::SweepMode::SingleSite;
		const double eSingle = single.CalculateFinite(sites, sweeps);

		DMRG::MPS::MPSDMRGAlgorithm two(2, ToTensor(Sz12), ToTensor(Sp12), 1., 1., m);
		two.sweepMode = DMRG::MPS::SweepMode::TwoSite;
		const double eTwo = two.CalculateFinite(sites, sweeps);

		std::cout << "  E(single-site) = " << eSingle << "   E(two-site) = " << eTwo << '\n';
		Check(cr, std::isfinite(eTwo), "two-site energy is finite");
		Check(cr, eTwo <= eSingle + 1e-6, "two-site energy <= single-site energy");
		Check(cr, static_cast<int>(two.results.size()) == sites - 1, "two-site chart has N-1 bond points");
	}
	std::cout << '\n';

	// ---- Single-site with density-matrix (subspace) perturbation --------------
	// The subspace-expansion correction should converge to essentially the same
	// energy as plain single-site (it helps escape local minima, not lower the
	// true variational minimum).
	std::cout << "-- Single-site subspace expansion (spin-1/2, 20 sites) --\n";
	{
		const int sites = 20, sweeps = 5;
		const unsigned int m = 20;

		DMRG::MPS::MPSDMRGAlgorithm plain(2, ToTensor(Sz12), ToTensor(Sp12), 1., 1., m);
		plain.sweepMode = DMRG::MPS::SweepMode::SingleSite;
		const double ePlain = plain.CalculateFinite(sites, sweeps);

		DMRG::MPS::MPSDMRGAlgorithm exp(2, ToTensor(Sz12), ToTensor(Sp12), 1., 1., m);
		exp.sweepMode = DMRG::MPS::SweepMode::SingleSiteSubspaceExpansion;
		exp.perturbationFactor = 1e-2;
		const double eExp = exp.CalculateFinite(sites, sweeps);

		std::cout << "  E(plain) = " << ePlain << "   E(subspace-expansion) = " << eExp << '\n';
		Check(cr, std::isfinite(eExp), "subspace-expansion energy is finite");
		Check(cr, std::abs(eExp - ePlain) < 5e-3, "subspace-expansion energy ~= single-site energy");
	}
	std::cout << '\n';

	// ---- Exact energy gap on the two-site spin-1/2 chain ----------------------
	// Eigenvalues of H = S1.S2 are: singlet -3/4 (ground) and triplet +1/4
	// (3-fold degenerate). The gap is therefore exactly 1.
	std::cout << "-- Two-site spin-1/2 energy gap (exact = 1.0) --\n";
	{
		for (int method = 0; method < 3; ++method)
		{
			std::cout << "  Method " << method << ": "
				<< (method == 0 ? "single-site" : (method == 1 ? "one-site with subspace expansion" : "two-site")) << '\n';
			DMRG::MPS::MPSDMRGAlgorithm dmrg(2, ToTensor(Sz12), ToTensor(Sp12), 1., 1., 16, method);
			dmrg.nrStates = 1;
			const double e0 = dmrg.CalculateFinite(2, 10);
			std::cout << "  E0 = " << e0 << "   gap = " << dmrg.EnergyGap << '\n';
			Check(cr, std::abs(e0 - (-0.75)) < 1e-4, "ground energy ~= -0.75");
			Check(cr, std::isfinite(dmrg.EnergyGap) && std::abs(dmrg.EnergyGap - 1.0) < 1e-3, "energy gap ~= 1.0");
		}
	}
	std::cout << '\n';

	// ---- Multiple excited states (loop over nrStates) -------------------------
	// The three triplet states are degenerate at +0.25, so every excited state
	// found on the two-site chain must have energy ~= +0.25.
	std::cout << "-- Multiple excited states (two-site spin-1/2) --\n";
	{
		for (int method = 0; method < 3; ++method)
		{
			std::cout << "  Method " << method << ": "
				<< (method == 0 ? "single-site" : (method == 1 ? "one-site with subspace expansion" : "two-site")) << '\n';
			DMRG::MPS::MPSDMRGAlgorithm dmrg(2, ToTensor(Sz12), ToTensor(Sp12), 1., 1., 16, method);
			dmrg.nrStates = 3;
			dmrg.CalculateFinite(2, 12);
			std::cout << "  excited energies:";
			for (double e : dmrg.excitedEnergies) std::cout << ' ' << e;
			std::cout << '\n';
			Check(cr, dmrg.excitedEnergies.size() == 3, "computed 3 excited states");
			bool allAtTriplet = dmrg.excitedEnergies.size() == 3;
			for (double e : dmrg.excitedEnergies)
				if (std::abs(e - 0.25) > 5e-3) allAtTriplet = false;
			Check(cr, allAtTriplet, "all excited energies ~= +0.25 (triplet)");
		}
	}
	std::cout << '\n';

	// ---- RNG seed independence ------------------------------------------------
	// The ground energy must not depend on the (reproducible) random seed used to
	// build the initial state.
	std::cout << "-- RNG seed independence (spin-1/2, 16 sites) --\n";
	{
		for (int method = 0; method < 3; ++method)
		{
			std::cout << "  Method " << method << ": "
				<< (method == 0 ? "single-site" : (method == 1 ? "one-site with subspace expansion" : "two-site")) << '\n';
			DMRG::MPS::MPSDMRGAlgorithm a(2, ToTensor(Sz12), ToTensor(Sp12), 1., 1., 16, method);
			a.rngSeed = 42;
			const double ea = a.CalculateFinite(16, 6);

			DMRG::MPS::MPSDMRGAlgorithm b(2, ToTensor(Sz12), ToTensor(Sp12), 1., 1., 16);
			b.rngSeed = 12345;
			const double eb = b.CalculateFinite(16, 6);

			std::cout << "  E(seed 42) = " << ea << "   E(seed 12345) = " << eb << '\n';
			Check(cr, std::abs(ea - eb) < 1e-4, "ground energy independent of RNG seed");
		}
	}
	std::cout << '\n';

	std::cout << "-- Multiple excited states (60-site spin-1) --\n";

	for (int method = 0; method < 3; ++method)
	{
		std::cout << "  Method " << method << ": "
			<< (method == 0 ? "single-site" : (method == 1 ? "one-site with subspace expansion" : "two-site")) << '\n';

		auto time1 = std::chrono::high_resolution_clock::now();
		DMRG::MPS::MPSDMRGAlgorithm dmrg(3, ToTensor(Sz1), ToTensor(Sp1), 1., 1., 32, 4, method);

		dmrg.CalculateFinite(60, 8);
		std::cout << "  excited energies:";
		for (double e : dmrg.excitedEnergies) std::cout << ' ' << e;
		std::cout << '\n';
		Check(cr, dmrg.excitedEnergies.size() == 4, "computed 4 excited states");
		std::cout << "  Haldane gap = " << dmrg.EnergyGap << '\n';
		Check(cr, std::isfinite(dmrg.EnergyGap) && std::abs(dmrg.EnergyGap - 0.4105) < 0.05, "Haldane gap ~= 0.4105");
		//std::cout << "  ground energy = " << dmrg.groundEnergy << '\n';

		double firstEnergy = dmrg.excitedEnergies.size() > 0 ? dmrg.excitedEnergies[0] : 0.;
		double secondEnergy = dmrg.excitedEnergies.size() > 1 ? dmrg.excitedEnergies[1] : 0.;
		double thirdEnergy = dmrg.excitedEnergies.size() > 2 ? dmrg.excitedEnergies[2] : 0.;
		bool isTriplet = std::abs(thirdEnergy - secondEnergy) < 1e-3 && std::abs(thirdEnergy - firstEnergy) < 1e-3;
		Check(cr, isTriplet, "three excited states are degenerate (triplet)");

		auto time2 = std::chrono::high_resolution_clock::now();
		std::cout << '\n';

		std::cout << "  Total runtime for this last test: "
			<< std::chrono::duration_cast<std::chrono::milliseconds>(time2 - time1).count() / 1000.0
			<< " seconds\n" << std::endl;
	}

	std::cout << "==================================================\n";
	if (cr.passed)
		std::cout << " ALL CHECKS PASSED\n";
	else
		std::cout << " " << cr.failures << " CHECK(S) FAILED\n";
	std::cout << "==================================================\n";

	return cr.passed ? 0 : 1;
}
