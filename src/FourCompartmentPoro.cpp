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

namespace mpet
{
FourCompartmentPoro::
FourCompartmentPoro() {}

void
FourCompartmentPoro::
setMaterialConstants()
{
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

    // extracellular/CSF constants
    A_e     = 1.0;
    B_e     = 0.99;
    M_e     = B_e * K / A_e / (1 - B_e * A_e); // N/m^2
    k_e     = 1.4e-14;    // m^2
    mu_e    = 8.9e-4;     // Ns/m^2
    kappa_e = k_e / mu_e;   // m^4/Ns

    // flow constants
    Q_p     = 5.8e-9;     // m^3/s
    Q_o     = Q_p * 1.;   // m^3/s
    R       = 8.5e13;     // m^-3
    p_bp    = 650.0;      // N/m^2
    p_bpA   = 13.3e3;     // N/m^2 arterial blood pressure (100mmHg)
}

void
FourCompartmentPoro::
setGridConstants(
    int m
    )
{
    int i;
    // grid properties
    J       = m;        // number of grid points
    dr      = L / (J - 1);  // grid spacing
    r       = Eigen::VectorXd(J); // radius vector
    for ( i = 0; i < J; i++ ) {
        r[i] = rV + i * dr;
    }
}

void
FourCompartmentPoro::
setTimeConstants(
    double initial_time,
    double final_time,
    double dtSecs
    )
{
    dt      = dtSecs;                  // time step size
    tf      = final_time;          // final solution time
    t0      = initial_time;           // initial time
    t       = t0;                       // current time
    N       = int((tf - t0 + dt / 10.) / dt) + 1; // number of time steps
}

void
FourCompartmentPoro::
initializeSystem()
{
    int i;
    int j;
    int n = 5;
    // solution properties
    x       = Eigen::VectorXd(n * J);       // LHS (solution) vector
    b       = Eigen::VectorXd(n * J);       // RHS vector
    for ( i = 0; i < n * J; i++ ) {
        x[i]      = 0.0;
        b[i]      = 0.0;
    }
    A       = Eigen::MatrixXd(n * J,n * J);   // Derivative matrix
    for ( i = 0; i < n * J; i++ ) {
        for ( j = 0; j < n * J; j++ ) {
            A(i, j) = 0.0;
        }
    }
}

void
FourCompartmentPoro::
createTransientDataFile ()
{
    transientFileName = baseName + "_transient.dat";
    std::ofstream fs;
    fs.open( transientFileName.c_str(), std::fstream::out );
    fs << "T, U, Pa, Pc, Pe, Pv" << std::endl;
    fs.close();
}

void
FourCompartmentPoro::
buildSystem ()
{
    // == build the A matrix ==== //
    int i;
    double sDot_ac, sDot_ce, sDot_cv, sDot_ev;
    for ( i = 1; i < J - 1; i++ ) {
        // set transfer fluxes
        sDot_ac = gamma_ac * abs(x[i + J] - x[i + 2 * J]);
        sDot_ce = gamma_ce * abs(x[i + 2 * J] - x[i + 3 * J]);
        sDot_cv = gamma_cv * abs(x[i + 2 * J] - x[i + 4 * J]);
        sDot_ev = gamma_ev * abs(x[i + 3 * J] - x[i + 4 * J]);

        // ** Displacement equation
        // section 1
        A(i, i - 1)        = 1. / pow(dr,2) - 1. / r[i] / dr;
        A(i, i)          = -2. / pow(dr,2) - 2. / pow(r[i],2);
        A(i, i + 1)        = 1. / pow(dr,2) + 1. / r[i] / dr;
        // section 2
        A(i, J + i - 1)      = (1. - 2. * nu) * A_a / 4. / G / (1 - nu) / dr;
        A(i, J + i + 1)      = -(1. - 2. * nu) * A_a / 4. / G / (1 - nu) / dr;
        // section 3
        A(i, 2 * J + i - 1)    = (1. - 2. * nu) * A_c / 4. / G / (1 - nu) / dr;
        A(i, 2 * J + i + 1)    = -(1. - 2. * nu) * A_c / 4. / G / (1 - nu) / dr;
        // section 4
        A(i, 3 * J + i - 1)    = (1. - 2. * nu) * A_e / 4. / G / (1 - nu) / dr;
        A(i, 3 * J + i + 1)    = -(1. - 2. * nu) * A_e / 4. / G / (1 - nu) / dr;
        // section 5
        A(i, 4 * J + i - 1)    = (1. - 2. * nu) * A_v / 4. / G / (1 - nu) / dr;
        A(i, 4 * J + i + 1)    = -(1. - 2. * nu) * A_v / 4. / G / (1 - nu) / dr;

        // ** arteriol pressure equation
        // section 5
        // section 6
        A(J + i, J + i - 1)    = 1. / pow(dr,2) - 1. / r[i] / dr;
        A(J + i, J + i)      = -2. / pow(dr,2);
        A(J + i, J + i + 1)    = 1. / pow(dr,2) + 1. / r[i] / dr;
        // section 7
        // section 8
        // b vector
        b[J + i]           = 1. / kappa_a * (sDot_ac);

        // ** capillary pressure equation
        // section 5
        // section 6
        A(2 * J + i, 2 * J + i - 1) = 1. / pow(dr,2) - 1. / r[i] / dr;
        A(2 * J + i, 2 * J + i)  = -2. / pow(dr,2);
        A(2 * J + i, 2 * J + i + 1) = 1. / pow(dr,2) + 1. / r[i] / dr;
        // section 7
        // section 8
        // b vector
        b[2 * J + i]         = 1. / kappa_c * (-sDot_ac + sDot_ce + sDot_cv);

        // ** CSF pressure equation
        // section 9
        // section 10
        // section 11
        A(3 * J + i, 3 * J + i - 1) = 1. / pow(dr,2) - 1. / r[i] / dr;
        A(3 * J + i, 3 * J + i)  = -2. / pow(dr,2);
        A(3 * J + i, 3 * J + i + 1) = 1. / pow(dr,2) + 1. / r[i] / dr;
        // section 12
        // b vector
        b[3 * J + i]         = 1. / kappa_e * (-sDot_ce + sDot_ev);

        // ** venous pressure equation
        // section 13
        // section 14
        // section 15
        // section 16
        A(4 * J + i, 4 * J + i - 1) = 1. / pow(dr,2) - 1. / r[i] / dr;
        A(4 * J + i, 4 * J + i)  = -2. / pow(dr,2);
        A(4 * J + i, 4 * J + i + 1) = 1. / pow(dr,2) + 1. / r[i] / dr;
        // b vector
        b[4 * J + i]         = 1. / kappa_v * (-sDot_cv - sDot_ev);
    }

    // Boundary Conditions

    // no displacement at skull
    A(J - 1, J - 1)          = 1.;

    // constant arterial blood pressure at the skull
    A(2 * J - 1, 2 * J - 1)      = 1.;
    b[2 * J - 1]             = p_bpA;

    // no capillary flow at the skull
    A(3 * J - 1, 3 * J - 2)      = -1.;
    A(3 * J - 1, 3 * J - 1)      = 1.;

    // CSF pressure at skull
    A(4 * J - 1, 4 * J - 1)      = 1.;
    b[4 * J - 1]             = p_bp + mu_e * R * Q_o;

    // constant venous blood pressure at the skull
    A(5 * J - 1, 5 * J - 1)      = 1.;
    b[5 * J - 1]             = p_bp;

    // stress equilibrium in the ventricle wall
    A(0, 0)              = 2. * nu / r[0] - (1. - nu) / dr;
    A(0, 1)              = (1. - nu) / dr;
    A(0, J)              = (1. - A_a) * (1. + nu) * (1. - 2. * nu) / E;
    A(0, 2 * J)            = (1. - A_c) * (1. + nu) * (1. - 2. * nu) / E;
    A(0, 3 * J)            = (1. - A_e) * (1. + nu) * (1. - 2. * nu) / E;
    A(0, 4 * J)            = (1. - A_v) * (1. + nu) * (1. - 2. * nu) / E;

    // no arteriol blood flow into ventricles
    A(J, J)              = -1.;
    A(J, J + 1)            = 1.;

    // capillary blood flow into ventricles
    A(2 * J, 2 * J)          = -1.;
    A(2 * J, 2 * J + 1)        = 1.;
    b[2 * J]               = mu_a * dr * Q_p / k_ce;

    // conservation of mass in ventricle
    A(3 * J, 0)            = 4.* M_PI* pow(r[0] + x[0],2) / dt;
    A(3 * J, 3 *
      J)          = M_PI * pow(d,4) / 128. / mu_e / L + 4.* M_PI* kappa_e* pow(
        r[0] + x[0],
        2) / dr;
    A(3 * J, 3 * J + 1)        = -4.* M_PI* kappa_e* pow(r[0] + x[0],2) / dr;
    A(3 * J, 4 * J - 1)        = -M_PI* pow(d,4) / 128. / mu_e / L;
    b[3 * J]               = Q_p + 4.* M_PI* pow(r[0] + x[0],2) * x[0] / dt;

    // venous blood flow into ventricles
    A(4 * J, 4 * J)          = -1.;
    A(4 * J, 4 * J + 1)        = 1.;
}

void
FourCompartmentPoro::
saveWallData ()
{
    wallFileName = baseName + "_wall.dat";
    std::ofstream wallFile( wallFileName.c_str());
    wallFile.precision(6);
    wallFile.setf(std::ios::scientific, std::ios::floatfield);
    wallFile << "r, u, p_a, p_c, p_e, p_v" << std::endl;
    int i;
    for ( i = 0; i < J; i++ ) {
        wallFile << r[i] << ", ";
        wallFile << x[i] << ", ";
        wallFile << x[i + J] << ", ";
        wallFile << x[i + 2 * J] << ", ";
        wallFile << x[i + 3 * J]  << ", ";
        wallFile << x[i + 4 * J] << std::endl;
    }
    wallFile.close();
}

void
FourCompartmentPoro::
saveTransientData ()
{
    std::ofstream fs;
    fs.open( transientFileName.c_str(), std::fstream::app );
    fs.precision(10);
    fs << t / 24. / 60. / 60. << ", ";
    fs.precision(6);
    fs << x[0] << ", ";
    fs << x[J] << ", ";
    fs << x[2 * J] << ", ";
    fs << x[3 * J] <<  ", ";
    fs << x[4 * J] << std::endl;
    fs.close();
}

void
FourCompartmentPoro::
Initialize (
    int m,
    double initialTime,
    double finalTime,
    double dtSecs,
    bool saveT,
    bool debugPrint,
    std::string bName
    )
{
    saveTrans = saveT;
    baseName = bName;
    mDebugPrint = debugPrint;
    // remove trailing spaces from the base name
    size_t endpos = baseName.find_last_not_of(" ");
    if ( std::string::npos != endpos ) {
        baseName = baseName.substr( 0, endpos + 1 );
    }

    setMaterialConstants();
    setGridConstants(m);
    setTimeConstants(initialTime, finalTime, dtSecs);
    initializeSystem();
    if (saveTrans) {
        createTransientDataFile();
    }
}

void
FourCompartmentPoro::
SetArteriolConstants(
    double Aa,
    double Ba,
    double kA,
    double muA
    )
{
    A_a     = Aa;
    B_a     = Ba;
    M_a     = B_a * K / A_a / (1 - B_a * A_a);
    k_a     = kA;
    mu_a    = muA;
    kappa_a = k_a / mu_a;
}

void
FourCompartmentPoro::
SetCapillaryConstants(
    double Ac,
    double Bc,
    double kC,
    double muC,
    double kCE
    )
{
    A_c     = Ac;
    B_c     = Bc;
    M_c     = B_c * K / A_c / (1 - B_c * A_c);
    k_c     = kC;
    mu_c    = muC;
    kappa_c = k_c / mu_c;
    k_ce    = kCE;
}

void
FourCompartmentPoro::
SetVenousConstants(
    double Av,
    double Bv,
    double kV,
    double muV
    )
{
    A_v     = Av;
    B_v     = Bv;
    M_v     = B_v * K / A_v / (1 - B_v * A_v);
    k_v     = kV;
    mu_v    = muV;
    kappa_v = k_v / mu_v;
}

void
FourCompartmentPoro::
SetTransferConstants(
    double gammaAC,
    double gammaCE,
    double gammaCV,
    double gammaEV
    )
{
    gamma_ac = gammaAC;
    gamma_ce = gammaCE;
    gamma_cv = gammaCV;
    gamma_ev = gammaEV;
}

void
FourCompartmentPoro::
Solve ()
{
    bool saveToFile = true;
    for (int i = 0; i < N; ++i) {
        t = t0 + i * dt;
        buildSystem();

        x = A.lu().solve(b);

        if (mDebugPrint) {
            double relative_error = (A*x - b).norm() / b.norm(); // norm() is L2 norm
            std::cout << "Current time: " << t << std::endl;
            std::cout << "... the relative error is:\n" << relative_error << std::endl;
        }

        if (saveTrans) {
            saveTransientData();
        }

        if (x[0] > 1e5) { // simulation is getting too big... cancel solution
            saveToFile = false;
            break;
        }
    }
    if (saveToFile) {
        saveWallData();
    }
}
}
