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
class FourCompartmentPoroOptions {
public:
    // arteriol constants
    double Aa;
    double Ba;
    double kA;
    double muA;
    // capillary constants
    double Ac;
    double Bc;
    double kC;
    double muC;
    double kCE;
    // venous constants
    double Av;
    double Bv;
    double kV;
    double muV;
    // transfer constants
    double gammaAC;
    double gammaCE;
    double gammaCV;
    double gammaEV;

    /**
     * Write the object to a stream
     *
     * @param os The output stream
     * @param rhs The vector to be written
     * @return
     */
    friend std::ostream& operator <<(
        std::ostream& os,
        const FourCompartmentPoroOptions& rhs
        )
    {
        return os << "FourCompartmentPoroOptions"
                  << std::endl << " Arteriol constants:"
                  << std::endl << "  - alpha: " << rhs.Aa
                  << std::endl << "  -  beta: " << rhs.Ba
                  << std::endl << "  - kappa: " << rhs.kA
                  << std::endl << "  -    mu: " << rhs.muA
                  << std::endl << " Capillary constants:"
                  << std::endl << "  - alpha: " << rhs.Ac
                  << std::endl << "  -  beta: " << rhs.Bc
                  << std::endl << "  - kappa: " << rhs.kC
                  << std::endl << "  -    mu: " << rhs.muC
                  << std::endl << "  -  k_ce: " << rhs.kCE
                  << std::endl << " Venous constants:"
                  << std::endl << "  - alpha: " << rhs.Av
                  << std::endl << "  -  beta: " << rhs.Bv
                  << std::endl << "  - kappa: " << rhs.kV
                  << std::endl << "  -    mu: " << rhs.muV
                  << std::endl << " Compartment transfer constants:"
                  << std::endl << "  - gamma_ac: " << rhs.gammaAC
                  << std::endl << "  - gamma_ce: " << rhs.gammaCE
                  << std::endl << "  - gamma_cv: " << rhs.gammaCV
                  << std::endl << "  - gamma_ev: " << rhs.gammaEV;
    }
};

class FourCompartmentPoro
{
public:
    // ************ public methods ************** //
    FourCompartmentPoro (
            int grid_size,
            double initial_time,
            double final_time,
            double dtSecs,
            bool saveTransientToFile,
            bool saveWallToFile,
            bool debugPrint,
            std::string bName,
            const FourCompartmentPoroOptions& opts);

    ~FourCompartmentPoro() {}

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
    bool mSaveTransientToFile;
    bool mSaveWallToFile;

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
