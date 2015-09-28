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
import csv
import multiprocessing as mp


class MonteCarloWorker(object):
    _STOP_SIGNAL = "STOP"

    def __init__(self, num_simulations=10000, blocked=False):
        self._num_simulations = num_simulations
        self._grid_size = 1500  # from resolution study
        self._num_steps = 3  # from resolution study
        secs_in_day = 86400.0
        self._initial_time = 0.0
        self._dt = secs_in_day / self._num_steps
        self._write_transient = False
        self._write_wall = False
        self._debug_print = False

        self._aqueduct_diameter = 0.25e-3 if blocked else 4e-3

        self._base_name = "monte_carlo"
        if blocked:
            self._base_name += "_blocked"
        self._out_filename = self._base_name + ".dat"
        self._out_csvfile = None

        self._num_procs = None
        self._num_simulations_per_proc = None
        self._run_q = None
        self._out_q = None
        self._out_p = None
        self._run_p = None

    def create_data(self, num_procs=None):
        """
        Run a series of simulations. As this is embarrasingly parallel
        we can use python multiprocess to compute them in parallel. The
        results of a run are pushed to the output queue to be written to
        file
        """
        self._num_procs = mp.cpu_count() if num_procs is None else num_procs
        self._num_simulations_per_proc = self._num_simulations / self._num_procs + 1
        self._run_q = mp.Queue()
        self._out_q = mp.Queue()
        self._run_p = [mp.Process(target=self._run_all, args=()) for i in range(self._num_procs)]
        self._out_p = mp.Process(target=self._output, args=())
        for p in self._run_p:
            p.start()
        self._out_p.start()
        for p in self._run_p:
            p.join()
        self._out_p.join()

    def load_results(self):
        """
        Open the output file and convert it to a numpy array ready for
        analysis
        """
        return np.loadtxt(self._out_filename, delimiter=',')

    def _run_all(self):
        """
        Runs all of the simulations, creating a unique id for each one
        """
        cur_proc = mp.current_process()._identity[0] - 1
        for i in range(self._num_simulations_per_proc):
            sim_idx = cur_proc * self._num_simulations_per_proc + i
            result = self._run(sim_idx)
            self._out_q.put(result)

        for i in range(self._num_procs):
            self._out_q.put(self._STOP_SIGNAL)

    def _output(self):
        """
        Take the results from each run and save it to file
        """
        with open(self._out_filename, "w") as outfile:
            self._out_csvfile = csv.writer(outfile)
            # Keep running until we see the stop message
            for works in range(self._num_procs):
                for result in iter(self._out_q.get, self._STOP_SIGNAL):
                    self._out_csvfile.writerow(result)

    def _run(self, sim_number):
        base_name = "{}_{:06d}".format(self._base_name, sim_number + 1)
        print base_name

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
        opts.aqueduct_diameter = self._aqueduct_diameter

        s = MPETSolver(self._grid_size,
                       self._initial_time,
                       self._num_steps,
                       self._dt,
                       self._write_transient,
                       self._write_wall,
                       self._debug_print,
                       base_name,
                       opts)
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


def main():
    from mpl_toolkits.mplot3d import Axes3D  # needed to register projection
    import matplotlib.pyplot as plt

    num_simulations = int(1e6)
    worker = MonteCarloWorker(num_simulations=num_simulations, blocked=False)
    rerun_simulations = True
    if rerun_simulations:
        num_procs = 3
        worker.create_data(num_procs=num_procs)
    results = worker.load_results()

    num_params = 11
    num_results = 5
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


if __name__ == "__main__":
    main()
