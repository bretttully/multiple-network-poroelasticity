"""!
@file
@date 07 Jun 2015

@license
Copyright 2015 Brett Tully

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""
import subprocess
import os
from mpet import FourCompartmentPoro
from mpet import profiling
from mpet import io


def main(run_profiling=False):
    grid_spacing = 81
    secs_in_day = 86400
    initial_time = 0.0
    final_time = 1.0 * secs_in_day
    dt = 100.0
    write_transient = False
    base_name = "example"
    debug_print = False

    # arteriol constants
    alpha_a = 1.0
    beta_a = 0.99
    k_a = 1e-10
    mu_a = 8.9e-4 * 3.  # about 3 times that of water

    # capillary constants
    alpha_c = 0.8
    beta_c = 0.99
    k_c = 1e-10
    mu_c = 8.9e-4 * 3.  # about 3 times that of water
    k_ce = 6e-4

    # venous constants
    alpha_v = 1.0
    beta_v = 0.99
    k_v = 1e-10
    mu_v = 8.9e-4 * 3.  # about 3 times that of water

    # transfer coefficients
    gamma_ac = 1.5e-19
    gamma_ce = 1.0e-22
    gamma_cv = 1.5e-19
    gamma_ev = 1.0e-13

    prof_prefix = os.path.join(os.getcwd(), base_name)
    if run_profiling:
        profiling.start(prof_prefix,
                        profile_memory=False,
                        cpu_profile_freq=1000)

    s = FourCompartmentPoro()
    s.initialize(grid_spacing, initial_time, final_time, dt, write_transient, debug_print, base_name)
    s.setArteriolConstants(alpha_a, beta_a, k_a, mu_a)
    s.setCapillaryConstants(alpha_c, beta_c, k_c, mu_c, k_ce)
    s.setVenousConstants(alpha_v, beta_v, k_v, mu_v)
    s.setTransferConstants(gamma_ac, gamma_ce, gamma_cv, gamma_ev)
    s.solve()

    if run_profiling:
        mem_profiling_on = profiling.profiling_memory_on()
        profiling.stop()
        prof_file = prof_prefix + ".prof"
        prof_pdf_file = os.path.join(prof_file + ".pdf")
        with open(prof_pdf_file, "w") as f:
            try:
                text = subprocess.check_output(["google-pprof",
                                                "--pdf",
                                                "/usr/bin/python",
                                                prof_file],
                                               stderr=subprocess.PIPE)
                f.write(text)
            except subprocess.CalledProcessError as err:
                print err
        if mem_profiling_on:
            heap_file = prof_prefix + ".0001.heap"
            heap_pdf_file = os.path.join(heap_file + ".pdf")
            with open(heap_pdf_file, "w") as f:
                try:
                    text = subprocess.check_output(["google-pprof",
                                                    "--pdf",
                                                    "--focus=hytrac",
                                                    "/usr/bin/python",
                                                    heap_file],
                                                   stderr=subprocess.PIPE)
                    f.write(text)
                except subprocess.CalledProcessError as err:
                    print err


if __name__ == "__main__":
    main()
