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
class FourCompartmentPoro
{
public:
    // ************ public methods ************** //
    FourCompartmentPoro ();
    ~FourCompartmentPoro() {}

    void initialize (
        int J,
        double initialTime,
        double finalTime,
        double dtSecs,
        bool saveTrans,
        bool debugPrint,
        std::string baseName
        );

    void setArteriolConstants(
        double A_a,
        double B_a,
        double k_a,
        double mu_a
        );

    void setCapillaryConstants(
        double A_c,
        double B_c,
        double k_c,
        double mu_c,
        double k_ce
        );

    void setVenousConstants(
        double A_v,
        double B_v,
        double k_v,
        double mu_v
        );

    void setTransferConstants(
        double gamma_ac,
        double gamma_ce,
        double gamma_cv,
        double gamma_ev
        );

    void solve();

private:
    // ************ private contants ************* //
    // geometry constants
    double rV;
    double rS;
    double L;
    double d;

    // poroelastic constants
    double nu;
    double E;
    double G;
    double K;

    // arteriol constants
    double A_a;
    double B_a;
    double M_a;
    double kappa_a;
    double k_a;
    double mu_a;

    // capillary constants
    double A_c;
    double B_c;
    double M_c;
    double kappa_c;
    double k_c;
    double mu_c;
    double k_ce;

    // extracellular/CSF constants
    double A_e;
    double B_e;
    double M_e;
    double kappa_e;
    double k_e;
    double mu_e;

    // venous constants
    double A_v;
    double B_v;
    double M_v;
    double kappa_v;
    double k_v;
    double mu_v;

    // transfer coefficients
    double gamma_ac;
    double gamma_ce;
    double gamma_cv;
    double gamma_ev;

    // flow contants
    double Q_p;
    double Q_o;
    double R;
    double p_bp;
    double p_bpA;

    // grid properties
    int J;
    Eigen::VectorXd r;
    double dr;

    // time control
    double dt;
    double tf;
    double t0;
    double t;
    int N;

    // solution properties
    static const int numElements = 5; // U, Pa, Pc, Pe, Pv
    Eigen::VectorXd x;
    Eigen::VectorXd b;
    Eigen::MatrixXd A;

    // output file constants
    std::string baseName;
    std::string wallFileName;
    std::string transientFileName;
    bool saveTrans;

    bool mDebugPrint;

    // ************ private methods ************** //
    void setMaterialConstants();
    void setGridConstants(
        int m
        );
    void setTimeConstants(
        double initialTime,
        double finalTime,
        double dtSecs
        );
    void initializeSystem();
    void buildSystem();
    void saveWallData() const;
    void createTransientDataFile() const;
    void saveTransientData() const;
};
}
