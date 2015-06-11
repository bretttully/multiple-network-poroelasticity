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

namespace
{
/**
 * Use this function to protect the wrapped objects from having additional
 * attributes assigned to them in the Python code.
 *
 *     .def("__setattr__", &protectSetAttr)
 *
 */
static void protectSetAttr(
    boost::python::object self,
    boost::python::str name,
    boost::python::object value
    )
{
    // if we can get the attribute, then we can set it as this means that
    // it has been set up by boost::python wrapping. If it hasn't been set
    // up, then this getattr attempt will throw an AttributeError
    boost::python::getattr(self, name);

    // now set the attribute using the Python C API. If you were to use
    // boost::python::setattr here, you would end up in an infinite loop
    PyObject_GenericSetAttr(self.ptr(), name.ptr(), value.ptr());
}
}

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
        .def("__setattr__", &protectSetAttr)
        ;

    class_<mpet::FourCompartmentPoroResult>("FourCompartmentPoroResult", init<>())
        .def(self_ns::str(self_ns::self))
        .def(self_ns::repr(self_ns::self))
        .def_readwrite("displacement", &mpet::FourCompartmentPoroResult::displacement)
        .def_readwrite("pressure_art", &mpet::FourCompartmentPoroResult::pressureArt)
        .def_readwrite("pressure_cap", &mpet::FourCompartmentPoroResult::pressureCap)
        .def_readwrite("pressure_csf", &mpet::FourCompartmentPoroResult::pressureCSF)
        .def_readwrite("pressure_ven", &mpet::FourCompartmentPoroResult::pressureVen)
        .def("__setattr__", &protectSetAttr)
        ;

    class_<mpet::FourCompartmentPoro>("FourCompartmentPoro",
                                      init<int, double, int, double,
                                           bool, bool, bool, std::string,
                                           const mpet::FourCompartmentPoroOptions&>())
        .def("solve", &mpet::FourCompartmentPoro::solve)
        .def("__setattr__", &protectSetAttr)
        ;
}
