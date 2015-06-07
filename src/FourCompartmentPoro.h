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
#include <string>
#include <iostream>
#include <fstream>
#include <math.h>
#include <eigen3/Eigen/Dense>

namespace mpet
{
class FourCompartmentPoro {

public:
    // ************ public methods ************** //
    FourCompartmentPoro ();
    ~FourCompartmentPoro() {}
    void Initialize (int J, double INIT_TIME_SECS, double FINAL_TIME_SECS, double DT_SECS, bool saveTrans, std::string baseName);    void SetArteriolConstants( double A_a, double B_a, double k_a, double mu_a );
    void SetCapillaryConstants( double A_c, double B_c, double k_c, double mu_c, double k_ce );
    void SetVenousConstants( double A_v, double B_v, double k_v, double mu_v );
    void SetTransferConstants( double gamma_ac, double gamma_ce, double gamma_cv, double gamma_ev );
    void Solve ();

private:
    // ************ private contants ************* //
    // geometry constants
    double rV, rS, L, d;

    // poroelastic constants
    double nu, E, G, K;

    // arteriol constants
    double A_a, B_a, M_a, kappa_a, k_a, mu_a;

    // capillary constants
    double A_c, B_c, M_c, kappa_c, k_c, mu_c, k_ce;

    // extracellular/CSF constants
    double A_e, B_e, M_e, kappa_e, k_e, mu_e;

    // venous constants
    double A_v, B_v, M_v, kappa_v, k_v, mu_v;

    // transfer coefficients
    double gamma_ac, gamma_ce, gamma_cv, gamma_ev;

    // flow contants
    double Q_p, Q_o, R, p_bp, p_bpA;

    // grid properties
    int J;
    Eigen::VectorXd r;
    double dr;

    // time control
    double dt, tf, t0, t;
    int N;

    // solution properties
    Eigen::VectorXd x;
    Eigen::VectorXd b;
    Eigen::MatrixXd A;

    // output file constants
    std::string baseName, transientFileName, wallFileName;
    bool saveTrans;

    // ************ private methods ************** //
    void setMaterialConstants ();
    void setGridConstants (int m);
    void setTimeConstants (double INIT_TIME_SECS, double FINAL_TIME_SECS, double DT_SECS);
    void initializeSystem ();
    void createTransientDataFile ();
    void buildSystem ();
    void saveWallData ();
    void saveTransientData ();
};
}
