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
from mpet import FourCompartmentPoro as CPPSolver
from mpet.four_compartment import FourCompartmentMPET as PythonSolver
from mpet import FourCompartmentPoroOptions
from mpet import profiling


def run(run_pure_python=False, run_profiling=False):
    grid_spacing = 81
    secs_in_day = 86400
    initial_time = 0.0
    final_time = 1.0 * secs_in_day
    num_steps = 1
    dt = (final_time - initial_time) / num_steps
    write_transient = False
    write_wall = True
    base_name = "example"
    if run_pure_python:
        base_name += "_python"
    else:
        base_name += "_cpp"
    debug_print = False

    opts = FourCompartmentPoroOptions()
    # arteriol constants
    opts.alpha_a = 1.0
    opts.beta_a = 0.99
    opts.kappa_a = 1e-10
    opts.mu_a = 8.9e-4 * 3.  # about 3 times that of water

    # capillary constants
    opts.alpha_c = 0.8
    opts.beta_c = 0.99
    opts.kappa_c = 1e-10
    opts.mu_c = 8.9e-4 * 3.  # about 3 times that of water
    opts.k_ce = 6e-4

    # venous constants
    opts.alpha_v = 1.0
    opts.beta_v = 0.99
    opts.kappa_v = 1e-10
    opts.mu_v = 8.9e-4 * 3.  # about 3 times that of water

    # transfer coefficients
    opts.gamma_ac = 1.5e-19
    opts.gamma_ce = 1.0e-22
    opts.gamma_cv = 1.5e-19
    opts.gamma_ev = 1.0e-13

    prof_prefix = os.path.join(os.getcwd(), base_name)
    if run_profiling:
        profiling.start(prof_prefix,
                        profile_memory=False,
                        cpu_profile_freq=1000)

    if run_pure_python:
        s = PythonSolver(grid_spacing, initial_time, final_time, dt,
                         write_transient, write_wall, debug_print, base_name, opts)
    else:
        s = CPPSolver(grid_spacing, initial_time, final_time, dt,
                      write_transient, write_wall, debug_print, base_name, opts)
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
    # import timeit
    # print "Pure Python:"
    # print timeit.timeit("from simple import run; run(run_pure_python=True)", number=10)
    # print "C++:"
    # print timeit.timeit("from simple import run; run(run_pure_python=False)", number=10)

    import matplotlib.pyplot as plt
    import numpy as np
    if True:
        import multiprocessing as mp
        pool = mp.Pool(processes=mp.cpu_count())
        pool.map(run, [True, False])
    else:
        run(run_pure_python=True)
        run(run_pure_python=False)

    p = np.genfromtxt("example_python_wall.dat", skiprows=1, delimiter=", ")
    c = np.genfromtxt("example_cpp_wall.dat", skiprows=1, delimiter=", ")
    plots = list()
    plots.append("Displacement")
    plots.append("Pressure: Art")
    plots.append("Pressure: Cap")
    plots.append("Pressure: CSF")
    plots.append("Pressure: Ven")
    num_plots = len(plots)
    for plt_idx, title in enumerate(plots):
        plt.subplot(2, num_plots, plt_idx + 1)
        plt.title(title)
        plt.plot(p[:, 0], p[:, plt_idx + 1], 'bx', label="Python")
        plt.plot(c[:, 0], c[:, plt_idx + 1], 'r-', label="C++")
        if plt_idx == 0:
            plt.legend()
        plt.subplot(2, num_plots, plt_idx + 1 + num_plots)
        plt.plot(p[:, 0], p[:, plt_idx + 6], 'bo')
        plt.plot(c[:, 0], c[:, plt_idx + 6], 'ro')
    plt.show()
