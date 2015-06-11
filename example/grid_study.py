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
__use_python_solver = True
if __use_python_solver:
    from mpet.four_compartment import FourCompartmentMPET as MPETSolver
else:
    from mpet import FourCompartmentPoro as MPETSolver
from mpet import FourCompartmentPoroOptions
from mpet import FourCompartmentPoroResult


class Setup(object):
    def __init__(self, grid_size, num_steps):
        self.grid_size = grid_size
        self.num_steps = num_steps

    def __str__(self):
        return "grid_size = {}, n_steps = {}".format(self.grid_size,
                                                     self.num_steps)


def run(setup):
    assert isinstance(setup, Setup)
    print "Running:", setup
    grid_size = setup.grid_size
    secs_in_day = 86400.0
    initial_time = 0.0
    num_steps = setup.num_steps
    dt = secs_in_day / num_steps
    write_transient = False
    write_wall = False
    base_name = "grid_study"
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

    # aqueduct diameter
    opts.aqueduct_diameter = 0.25e-3# 4e-3

    s = MPETSolver(grid_size, initial_time, num_steps, dt, write_transient,
                   write_wall, debug_print, base_name, opts)
    result = s.solve()
    assert isinstance(result, FourCompartmentPoroResult)
    return [setup.grid_size,
            setup.num_steps,
            result.displacement,
            result.pressure_art,
            result.pressure_cap,
            result.pressure_csf,
            result.pressure_ven]


if __name__ == "__main__":
    from mpl_toolkits.mplot3d import Axes3D  # needed to register projection
    import matplotlib.pyplot as plt
    import numpy as np
    import multiprocessing as mp
    fname = "grid_study.dat"

    rerun_simulations = False
    if rerun_simulations:
        grid_size_list = range(500, 5000 + 1, 200)
        num_steps_list = range(3, 4)
        setup_list = list()
        for g in grid_size_list:
            for n in num_steps_list:
                setup_list.append(Setup(g, n))

        pool = mp.Pool(processes=2)#mp.cpu_count())
        results = np.array(pool.map(run, setup_list))
        np.savetxt(fname, results, delimiter='\t')
    else:
        results = np.loadtxt(fname, delimiter='\t')

    plots = list()
    plots.append("Displacement")
    plots.append("Pressure: Art")
    plots.append("Pressure: Cap")
    plots.append("Pressure: CSF")
    plots.append("Pressure: Ven")
    num_plots = len(plots)
    fig = plt.figure()
    for plt_idx, plt_title in enumerate(plots):
        # ax = fig.add_subplot(1, num_plots, plt_idx + 1, projection='3d')
        # ax.scatter(results[:, 0], results[:, 1], zs=results[:, plt_idx + 2])
        ax = fig.add_subplot(1, num_plots, plt_idx + 1)
        ax.plot(results[:, 0], results[:, plt_idx + 2], 'o')
        ax.grid()
        ax.set_title(plt_title)
    plt.show()
