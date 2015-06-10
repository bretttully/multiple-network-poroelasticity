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

#include <boost/python.hpp>
#include <FourCompartmentPoro.h>

// the name of the module must match the filename
BOOST_PYTHON_MODULE(libpympet) {
    using namespace boost::python;

    class_<mpet::FourCompartmentPoroOptions>("FourCompartmentPoroOptions", init<>())
        .def(self_ns::str(self_ns::self))
        .def(self_ns::repr(self_ns::self))
        // arteriol constants
        .def_readwrite("alpha_a", &mpet::FourCompartmentPoroOptions::Aa)
        .def_readwrite("beta_a", &mpet::FourCompartmentPoroOptions::Ba)
        .def_readwrite("kappa_a", &mpet::FourCompartmentPoroOptions::kA)
        .def_readwrite("mu_a", &mpet::FourCompartmentPoroOptions::muA)
        // capillary constants
        .def_readwrite("alpha_c", &mpet::FourCompartmentPoroOptions::Ac)
        .def_readwrite("beta_c", &mpet::FourCompartmentPoroOptions::Bc)
        .def_readwrite("kappa_c", &mpet::FourCompartmentPoroOptions::kC)
        .def_readwrite("mu_c", &mpet::FourCompartmentPoroOptions::muC)
        .def_readwrite("k_ce", &mpet::FourCompartmentPoroOptions::kCE)
        // venous constants
        .def_readwrite("alpha_v", &mpet::FourCompartmentPoroOptions::Av)
        .def_readwrite("beta_v", &mpet::FourCompartmentPoroOptions::Bv)
        .def_readwrite("kappa_v", &mpet::FourCompartmentPoroOptions::kV)
        .def_readwrite("mu_v", &mpet::FourCompartmentPoroOptions::muV)
        // transfer constants
        .def_readwrite("gamma_ac", &mpet::FourCompartmentPoroOptions::gammaAC)
        .def_readwrite("gamma_ce", &mpet::FourCompartmentPoroOptions::gammaCE)
        .def_readwrite("gamma_cv", &mpet::FourCompartmentPoroOptions::gammaCV)
        .def_readwrite("gamma_ev", &mpet::FourCompartmentPoroOptions::gammaEV)
        // geometric properties
        .def_readwrite("aqueduct_diameter", &mpet::FourCompartmentPoroOptions::aqueductDiameter)
//        .def_pickle(ns::controllerResultPickleSuite())
        ;

    class_<mpet::FourCompartmentPoro>("FourCompartmentPoro",
                                      init<int, double, int, double,
                                           bool, bool, bool, std::string,
                                           const mpet::FourCompartmentPoroOptions&>())
        .def("solve", &mpet::FourCompartmentPoro::solve)
        ;
}
