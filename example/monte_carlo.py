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
import numpy as np


def run(sim_number):
    grid_size = 1500  # from resolution study
    num_steps = 3  # from resolution study
    secs_in_day = 86400.0
    initial_time = 0.0
    dt = secs_in_day / num_steps
    write_transient = False
    write_wall = False
    base_name = "monte_carlo_{:06d}".format(sim_number + 1)
    print base_name
    debug_print = False

    opts = FourCompartmentPoroOptions()
    # arteriol constants
    opts.alpha_a = np.random.uniform(0.8, 1.0)
    opts.beta_a = 0.99
    opts.kappa_a = 10 ** np.random.uniform(-14.0, -8.0)
    opts.mu_a = 8.9e-4 * 3.  # about 3 times that of water

    # capillary constants
    opts.alpha_c = np.random.uniform(0.8, 1.0)
    opts.beta_c = 0.99
    opts.kappa_c = 10 ** np.random.uniform(-14.0, -8.0)
    opts.mu_c = 8.9e-4 * 3.  # about 3 times that of water
    opts.k_ce = 10 ** np.random.uniform(-10.0, 0.0)

    # venous constants
    opts.alpha_v = np.random.uniform(0.8, 1.0)
    opts.beta_v = 0.99
    opts.kappa_v = 10 ** np.random.uniform(-14.0, -8.0)
    opts.mu_v = 8.9e-4 * 3.  # about 3 times that of water

    # transfer coefficients
    opts.gamma_ac = 10 ** np.random.uniform(np.log10(3.0e-19), np.log10(1.0e-19))
    opts.gamma_ce = 10 ** np.random.uniform(-20, np.log10(2.5e-19))
    opts.gamma_cv = 10 ** np.random.uniform(-20, np.log10(2.5e-19))
    opts.gamma_ev = 10 ** np.random.uniform(-13.0, -8.0)

    # aqueduct diameter
    opts.aqueduct_diameter = 0.25e-3# 4e-3

    s = MPETSolver(grid_size, initial_time, num_steps, dt, write_transient,
                   write_wall, debug_print, base_name, opts)
    result = s.solve()
    assert isinstance(result, FourCompartmentPoroResult)
    return [opts.alpha_a,
            opts.alpha_c,
            opts.alpha_v,
            opts.kappa_a,
            opts.kappa_c,
            opts.kappa_v,
            opts.k_ce,
            opts.gamma_ac,
            opts.gamma_ce,
            opts.gamma_cv,
            opts.gamma_ev,
            result.displacement,
            result.pressure_art,
            result.pressure_cap,
            result.pressure_csf,
            result.pressure_ven]


if __name__ == "__main__":
    from mpl_toolkits.mplot3d import Axes3D  # needed to register projection
    import matplotlib.pyplot as plt
    import multiprocessing as mp
    fname = "monte_carlo.dat"

    num_simulations = 10000
    rerun_simulations = True
    if rerun_simulations:
        pool = mp.Pool(processes=mp.cpu_count())
        results = np.array(pool.map(run, range(num_simulations)))
        np.savetxt(fname, results, delimiter='\t')
    else:
        results = np.loadtxt(fname, delimiter='\t')

    num_params = 11
    num_results = 5
    assert results.shape[0] == num_simulations
    assert results.shape[1] == num_params + num_results

    param_names = list()
    param_names.append("alpha_a")
    param_names.append("alpha_c")
    param_names.append("alpha_v")
    param_names.append("kappa_a")
    param_names.append("kappa_c")
    param_names.append("kappa_v")
    param_names.append("k_ce")
    param_names.append("gamma_ac")
    param_names.append("gamma_ce")
    param_names.append("gamma_cv")
    param_names.append("gamma_ev")
    assert len(param_names) == num_params

    result_names = list()
    result_names.append("Displacement")
    result_names.append("Pressure: Art")
    result_names.append("Pressure: Cap")
    result_names.append("Pressure: CSF")
    result_names.append("Pressure: Ven")
    assert len(result_names) == num_results

    fig = plt.figure()
    for param_idx, param_title in enumerate(param_names):
        for result_idx, result_title in enumerate(result_names):
            plt_idx = result_idx * num_params+ param_idx + 1
            ax = fig.add_subplot(num_results, num_params, plt_idx)
            if param_title.startswith("alpha"):
                ax.plot(results[:, param_idx], results[:, result_idx + num_params], 'o')
            else:
                ax.semilogx(results[:, param_idx], results[:, result_idx + num_params], 'o')
            ax.grid()
            if result_idx == 0:
                ax.set_title(param_title)
            if param_idx == 0:
                ax.set_ylabel(result_title)
    plt.show()
