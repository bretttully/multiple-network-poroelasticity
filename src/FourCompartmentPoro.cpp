/**
 * @file
 * @date 07 Jun 2015
 *
 * @license
 * Copyright 2015 Brett Tully
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#include "FourCompartmentPoro.h"
#include <eigen3/Eigen/Sparse>

namespace mpet
{
FourCompartmentPoro::
FourCompartmentPoro(
    int grid_size,
    double initial_time,
    double final_time,
    double dtSecs,
    bool saveTransientToFile,
    bool saveWallToFile,
    bool debugPrint,
    std::string bName,
    const FourCompartmentPoroOptions& opts
    )
{
    mSaveTransientToFile = saveTransientToFile;
    mSaveWallToFile = saveWallToFile;
    baseName = bName;
    mDebugPrint = debugPrint;
    // remove trailing spaces from the base name
    size_t endpos = baseName.find_last_not_of(" ");
    if ( std::string::npos != endpos ) {
        baseName = baseName.substr( 0, endpos + 1 );
    }
    transientFileName = baseName + "_transient.dat";
    wallFileName = baseName + "_wall.dat";

    // geometry constants
    rV      = 30.0e-3;    // m
    rS      = 100.0e-3;   // m
    L       = rS - rV;    // m
    d       = 4e-3; //0.25e-3;    // m (assume a blocked aqueduct 0.25mm)

    // poroelastic constants
    E       = 584.0;      // N/m^2
    nu      = 0.35;
    G       = E / 2. / (1. + nu); // N/m^2
    K       = E / 3. / (1. - 2. * nu); // N/m^2

    // arteriol constants
    A_a     = opts.Aa;
    B_a     = opts.Ba;
    M_a     = B_a * K / A_a / (1 - B_a * A_a);
    k_a     = opts.kA;
    mu_a    = opts.muA;
    kappa_a = k_a / mu_a;

    // capillary constants
    A_c     = opts.Ac;
    B_c     = opts.Bc;
    M_c     = B_c * K / A_c / (1 - B_c * A_c);
    k_c     = opts.kC;
    mu_c    = opts.muC;
    kappa_c = k_c / mu_c;
    k_ce    = opts.kCE;

    // extracellular/CSF constants
    A_e     = 1.0;
    B_e     = 0.99;
    M_e     = B_e * K / A_e / (1 - B_e * A_e); // N/m^2
    k_e     = 1.4e-14;    // m^2
    mu_e    = 8.9e-4;     // Ns/m^2
    kappa_e = k_e / mu_e;   // m^4/Ns

    // venous constants
    A_v     = opts.Av;
    B_v     = opts.Bv;
    M_v     = B_v * K / A_v / (1 - B_v * A_v);
    k_v     = opts.kV;
    mu_v    = opts.muV;
    kappa_v = k_v / mu_v;

    // compartment transfer constants
    gamma_ac = opts.gammaAC;
    gamma_ce = opts.gammaCE;
    gamma_cv = opts.gammaCV;
    gamma_ev = opts.gammaEV;

    // flow constants
    Q_p     = 5.8e-9;     // m^3/s
    Q_o     = Q_p * 1.;   // m^3/s
    R       = 8.5e13;     // m^-3
    p_bp    = 650.0;      // N/m^2
    p_bpA   = 13.3e3;     // N/m^2 arterial blood pressure (100mmHg)

    // grid properties
    J       = grid_size;          // number of grid points
    dr      = L / (J - 1);        // grid spacing
    r       = Eigen::VectorXd(J); // radius vector
    for (int i = 0; i < J; ++i) {
        r[i] = rV + i * dr;
    }
    x = Eigen::VectorXd(numElements * J);                 // LHS (solution) vector
    b = Eigen::VectorXd(numElements * J);                 // RHS vector
    A = Eigen::MatrixXd(numElements * J,numElements * J); // Derivative matrix
    residual = Eigen::VectorXd(numElements * J);          // r = A*x - b
    x.setZero();
    b.setZero();
    A.setZero();
    residual.setZero();

    // simulation properties
    dt      = dtSecs;                             // time step size
    tf      = final_time;                         // final solution time
    t0      = initial_time;                       // initial time
    t       = t0;                                 // current time
    N       = int((tf - t0 + dt / 10.) / dt) + 1; // number of time steps
}

void
FourCompartmentPoro::
initialiseSystem()
{
    // == build the A matrix ==== //
    double oneOverDrSquared = 1.0 / (dr * dr);
    double strainMultiplier = (1.0 - 2.0 * nu) / (4.0 * G * (1.0 - nu) * dr);
    for (int i = 1; i < J - 1; ++i) {
        double displacementMultiplierMinus = oneOverDrSquared - 1. / (r[i] * dr);
        double displacementMultiplierCenter = -2. * oneOverDrSquared;
        double displacementMultiplierPlus = oneOverDrSquared + 1. / (r[i] * dr);

        // n = 0: Displacement equation
        // n = 1: arteriol pressure equation
        // n = 2: capillary pressure equation
        // n = 3: CSF pressure equation
        // n = 4: venous pressure equation
        for (int n = 0; n < 5; ++n) {
            A(n * J + i, n * J + i - 1) = displacementMultiplierMinus;
            A(n * J + i, n * J + i) = displacementMultiplierCenter;
            A(n * J + i, n * J + i + 1) = displacementMultiplierPlus;
        }

        // update displacement equation
        A(i, i) -= 2.0 / std::pow(r[i], 2);

        // section 2
        A(i, J + i - 1) = strainMultiplier * A_a;
        A(i, J + i + 1) = -strainMultiplier * A_a;
        // section 3
        A(i, 2 * J + i - 1) = strainMultiplier * A_c;
        A(i, 2 * J + i + 1) = -strainMultiplier * A_c;
        // section 4
        A(i, 3 * J + i - 1) = strainMultiplier * A_e;
        A(i, 3 * J + i + 1) = -strainMultiplier * A_e;
        // section 5
        A(i, 4 * J + i - 1) = strainMultiplier * A_v;
        A(i, 4 * J + i + 1) = -strainMultiplier * A_v;
    }

    // ---
    // Boundary Conditions

    // no displacement at skull
    A(J - 1, J - 1) = 1.0;

    // constant arterial blood pressure at the skull
    A(2 * J - 1, 2 * J - 1) = 1.0;
    b[2 * J - 1] = p_bpA;

    // no capillary flow at the skull
    A(3 * J - 1, 3 * J - 2) = -1.0;
    A(3 * J - 1, 3 * J - 1) = 1.0;

    // CSF pressure at skull
    A(4 * J - 1, 4 * J - 1) = 1.0;
    b[4 * J - 1] = p_bp + mu_e * R * Q_o;

    // constant venous blood pressure at the skull
    A(5 * J - 1, 5 * J - 1) = 1.0;
    b[5 * J - 1] = p_bp;

    // stress equilibrium in the ventricle wall
    A(0, 0) = 2.0 * nu / r[0] - (1.0 - nu) / dr;
    A(0, 1) = (1.0 - nu) / dr;
    double stressMultiplier = (1.0 + nu) * (1.0 - 2.0 * nu) / E;
    A(0, 1 * J) = (1. - A_a) * stressMultiplier;
    A(0, 2 * J) = (1. - A_c) * stressMultiplier;
    A(0, 3 * J) = (1. - A_e) * stressMultiplier;
    A(0, 4 * J) = (1. - A_v) * stressMultiplier;

    // no arteriol blood flow into ventricles
    A(J, J) = -1.0;
    A(J, J + 1) = 1.0;

    // capillary blood flow into ventricles
    A(2 * J, 2 * J) = -1.0;
    A(2 * J, 2 * J + 1) = 1.0;
    b[2 * J] = mu_a * dr * Q_p / k_ce;

    // venous blood flow into ventricles
    A(4 * J, 4 * J) = -1.0;
    A(4 * J, 4 * J + 1) = 1.0;
}

void
FourCompartmentPoro::
updateTimeDependentSystem()
{
    // == build the A matrix ==== //
    double sDot_ac, sDot_ce, sDot_cv, sDot_ev;
    for (int i = 1; i < J - 1; ++i) {
        // set transfer fluxes
        sDot_ac = gamma_ac * std::fabs(x[i + 1 * J] - x[i + 2 * J]);
        sDot_ce = gamma_ce * std::fabs(x[i + 2 * J] - x[i + 3 * J]);
        sDot_cv = gamma_cv * std::fabs(x[i + 2 * J] - x[i + 4 * J]);
        sDot_ev = gamma_ev * std::fabs(x[i + 3 * J] - x[i + 4 * J]);

        // arteriol pressure equation
        b[J + i] = (sDot_ac) / kappa_a;
        // capillary pressure equation
        b[2 * J + i] = (-sDot_ac + sDot_ce + sDot_cv) / kappa_c;
        // CSF pressure equation
        b[3 * J + i] = (-sDot_ce + sDot_ev) / kappa_e;
        // venous pressure equation
        b[4 * J + i] = (-sDot_cv - sDot_ev) / kappa_v;
    }

    // Boundary Conditions

    // conservation of mass in ventricle
    double d2 = d * d;
    double d4 = d2 * d2;
    double const_1 = M_PI * d4 / (128. * mu_e * L);
    double const_2 = 4.0 * M_PI * (r[0] + x[0]) * (r[0] + x[0]);
    A(3 * J, 0) = const_2 / dt;
    A(3 * J, 3 *J) = const_1 + const_2 * kappa_e / dr;
    A(3 * J, 3 * J + 1) = -const_2 * kappa_e / dr;
    A(3 * J, 4 * J - 1) = -const_1;
    b[3 * J] = Q_p + const_2 * x[0] / dt;
}

void
FourCompartmentPoro::
saveWallData() const
{
    std::ofstream fs(wallFileName.c_str());
    fs.precision(6);
    fs.setf(std::ios::scientific, std::ios::floatfield);
    fs << "r, u, p_a, p_c, p_e, p_v"
       << "u_res, p_a_res, p_c_res, p_e_res, p_v_res"
       << std::endl;
    int i;
    for ( i = 0; i < J; i++ ) {
        fs << r[i];
        fs << ", " << x[i + 0 * J];
        fs << ", " << x[i + 1 * J];
        fs << ", " << x[i + 2 * J];
        fs << ", " << x[i + 3 * J];
        fs << ", " << x[i + 4 * J];
        fs << ", " << residual[i + 0 * J];
        fs << ", " << residual[i + 1 * J];
        fs << ", " << residual[i + 2 * J];
        fs << ", " << residual[i + 3 * J];
        fs << ", " << residual[i + 4 * J];
        fs << std::endl;
    }
    fs.close();
}

void
FourCompartmentPoro::
createTransientDataFile() const
{
    std::ofstream fs;
    fs.open(transientFileName.c_str(), std::fstream::out );
    fs << "T, U, Pa, Pc, Pe, Pv, "
       << "U_res, Pa_res, Pc_res, Pe_res, Pv_res"
       << std::endl;
    fs.close();
}

void
FourCompartmentPoro::
saveTransientData() const
{
    std::ofstream fs;
    fs.open( transientFileName.c_str(), std::fstream::app );
    fs.precision(10);
    fs << t / 24. / 60. / 60.;
    fs.precision(6);
    fs << ", " << x[0 * J];
    fs << ", " << x[1 * J];
    fs << ", " << x[2 * J];
    fs << ", " << x[3 * J];
    fs << ", " << x[4 * J];
    fs << ", " << residual[0 * J];
    fs << ", " << residual[1 * J];
    fs << ", " << residual[2 * J];
    fs << ", " << residual[3 * J];
    fs << ", " << residual[4 * J];
    fs << std::endl;
    fs.close();
}

void
FourCompartmentPoro::
solve()
{
    if (mSaveTransientToFile) {
        createTransientDataFile();
    }

    initialiseSystem();

    bool simulationFailed = false;
    for (int i = 0; !simulationFailed && i < N; ++i) {
        t = t0 + i * dt;
        updateTimeDependentSystem();

//        using SpMat = Eigen::SparseMatrix<double>;
//        SpMat sparseA = A.sparseView();
//        Eigen::SimplicialCholesky<SpMat> chol(sparseA);  // performs a Cholesky factorization of A
//        x = chol.solve(b);         // use the factorization to solve for the given

        x = A.lu().solve(b);
//        x = A.householderQr().solve(b);
        residual = A * x - b;

        if (mDebugPrint) {
            double relative_error = residual.norm() / b.norm(); // norm() is L2 norm
            std::cout << "Current time: " << t << " sec. The relative error is: " << relative_error << std::endl;
        }

        if (mSaveTransientToFile) {
            saveTransientData();
        }

        if (x[0] > 1e5) { // simulation is getting too big... cancel solution
            simulationFailed = true;
        }
    }
    if (mSaveWallToFile && !simulationFailed) {
        saveWallData();
    }
}
}
