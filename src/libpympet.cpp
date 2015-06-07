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

    class_<mpet::FourCompartmentPoro>("FourCompartmentPoro", init<>())
//        .def(self_ns::str(self_ns::self))
//        .def(self_ns::repr(self_ns::self))
    .def("initialize", &mpet::FourCompartmentPoro::Initialize)
    .def("setArteriolConstants",
         &mpet::FourCompartmentPoro::SetArteriolConstants)
    .def("setCapillaryConstants",
         &mpet::FourCompartmentPoro::SetCapillaryConstants)
    .def("setVenousConstants", &mpet::FourCompartmentPoro::SetVenousConstants)
    .def("setTransferConstants",
         &mpet::FourCompartmentPoro::SetTransferConstants)
    .def("solve", &mpet::FourCompartmentPoro::Solve)
//        .def_pickle(ns::controllerResultPickleSuite())
    ;
}
