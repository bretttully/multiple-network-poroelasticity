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
import os
import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spsl
import matplotlib.pyplot as plt
from mpet import FourCompartmentPoroOptions


class FourCompartmentMPET(object):
    def __init__(self, grid_spacing, initial_time, num_steps, dt,
                 write_transient, write_wall, debug_print, base_name, opts):
        assert isinstance(opts, FourCompartmentPoroOptions)

        # geometry constants
        self.rV = 30.0e-3  # m
        self.rS = 100.0e-3  # m
        self.L = self.rS - self.rV  # m
        self.d = opts.aqueduct_diameter

        # poroelastic constants
        self.nu = 0.35
        self.E = 584.0  # N/m^2
        self.G = self.E / (2. * (1. + self.nu))  # N/m^2
        self.K = self.E / (3. * (1. - 2. * self.nu))  # N/m^2

        # arteriol constants
        self.A_a = opts.alpha_a
        self.B_a = opts.beta_a
        self.M_a = self.B_a * self.K / self.A_a / (1 - self.B_a * self.A_a)  # N/m^2
        self.k_a = opts.kappa_a  # m^2
        self.mu_a = opts.mu_a  # Ns/m^2
        self.kappa_a = self.k_a / self.mu_a  # m^4/Ns

        # capillary constants
        self.A_c = opts.alpha_c
        self.B_c = opts.beta_c
        self.M_c = self.B_c * self.K / self.A_c / (1 - self.B_c * self.A_c)  # N/m^2
        self.k_c = opts.kappa_c  # m^2
        self.mu_c = opts.mu_c  # Ns/m^2
        self.kappa_c = self.k_c / self.mu_c  # m^4/Ns
        self.k_ce = opts.k_ce

        # extracellular/CSF constants
        self.A_e = 1.0
        self.B_e = 0.99
        self.M_e = self.B_e * self.K / self.A_e / (1 - self.B_e * self.A_e)  # N/m^2
        self.k_e = 1.4e-14  # m^2
        self.mu_e = 8.9e-4  # Ns/m^2
        self.kappa_e = self.k_e / self.mu_e  # m^4/Ns

        # venous constants
        self.A_v = opts.alpha_v
        self.B_v = opts.beta_v
        self.M_v = self.B_v * self.K / self.A_v / (1 - self.B_v * self.A_v)  # N/m^2
        self.k_v = opts.kappa_v  # m^2
        self.mu_v = opts.mu_v  # Ns/m^2
        self.kappa_v = self.k_v / self.mu_v  # m^4/Ns

        # transfer coefficients
        self.gamma_ac = opts.gamma_ac
        self.gamma_ce = opts.gamma_ce
        self.gamma_cv = opts.gamma_cv
        self.gamma_ev = opts.gamma_ev

        # flow contants
        self.Q_p = 5.8e-9  # m^3/s
        self.Q_o = self.Q_p * 1.0  # m^3/s
        self.R = 8.5e13  # m^-3
        self.p_bp = 650.0  # N/m^2
        self.p_bpA = 13.3e3  # N/m^2 arterial blood pressure (100mmHg)

        # grid properties
        self.J = grid_spacing
        self.r = np.linspace(self.rV, self.rS, self.J)
        self.dr = self.L / (self.J - 1)
        n = 5  # U, Pa, Pc, Pe, Pv
        vec_size = n * self.J
        self.x = np.zeros((vec_size,))  # LHS (solution) vector
        self.b = None  # RHS vector
        self.A = None  # Derivative matrix
        self.residual = None  # residual = A*x - b

        # time control
        self.dt = dt  # time step size
        self.t0 = initial_time  # initiali solution time
        self.t = self.t0  # current time
        self.N = num_steps  # number of time steps

        # output file constants
        self.base_name = base_name
        self.transient_fileName = base_name + "_transient.dat"
        self.wall_file_name = base_name + "_wall.dat"
        self.save_transient = write_transient
        self.save_wall = write_wall
        self.debug_print = debug_print

    def _build_system(self):
        # == build the A & b matrix ==
        vec_size = self.x.shape[0]
        self.b = np.zeros((vec_size,))
        self.A = np.zeros((vec_size, vec_size))

        # ---
        # useful consts
        A = self.A
        b = self.b
        x = self.x
        J = self.J
        dr = self.dr
        one_over_dr_squared = 1.0 / (dr * dr)
        r = self.r
        central_r = self.r[1:-1]

        # ---
        # set up the difference tridiagonals
        difference_eye_minus = np.diag(one_over_dr_squared - 1. / (r[1:] * dr), -1)
        difference_eye_central = -2.0 * one_over_dr_squared * np.eye(J, k=0)
        difference_eye_plus = np.diag(one_over_dr_squared + 1. / (r[:-1] * dr), 1)
        difference_tridiagonal = difference_eye_minus + difference_eye_central + difference_eye_plus
        difference_tridiagonal[0, :] = 0.0
        difference_tridiagonal[-1, :] = 0.0
        for i in range(self.x.shape[0] / self.J):
            A[i * J:(i + 1) * J, i * J:(i + 1) * J] = difference_tridiagonal
        # update for the displacement equation
        A[1:J - 1, 1:J - 1] -= np.diag(2. / (central_r * central_r), 0)

        # ---
        # set transfer fluxes
        sDot_ac = self.gamma_ac * np.fabs(x[1 * J + 1:2 * J - 1] - x[2 * J + 1:3 * J - 1])
        sDot_ce = self.gamma_ce * np.fabs(x[2 * J + 1:3 * J - 1] - x[3 * J + 1:4 * J - 1])
        sDot_cv = self.gamma_cv * np.fabs(x[2 * J + 1:3 * J - 1] - x[4 * J + 1:5 * J - 1])
        sDot_ev = self.gamma_ev * np.fabs(x[3 * J + 1:4 * J - 1] - x[4 * J + 1:5 * J - 1])
        # arteriol pressure equation
        b[J + 1:2 * J - 1] = sDot_ac / self.kappa_a
        # capillary pressure equation
        b[2 * J + 1:3 * J - 1] = (-sDot_ac + sDot_ce + sDot_cv) / self.kappa_c
        # CSF pressure equation
        b[3 * J + 1:4 * J - 1] = (-sDot_ce + sDot_ev) / self.kappa_e
        # venous pressure equation
        b[4 * J + 1:5 * J - 1] = (-sDot_cv - sDot_ev) / self.kappa_v

        # ---
        # set up the poro diagonals
        strain_multiplier = (1.0 - 2.0 * self.nu) / (4.0 * self.G * (1.0 - self.nu) * dr)
        displacement_tridiagonal = (np.eye(J, k=-1) - np.eye(J, k=1)) * strain_multiplier
        displacement_tridiagonal[0, :] = 0.0
        displacement_tridiagonal[-1, :] = 0.0
        A[:J, 1 * J:2 * J] = displacement_tridiagonal * self.A_a
        A[:J, 2 * J:3 * J] = displacement_tridiagonal * self.A_c
        A[:J, 3 * J:4 * J] = displacement_tridiagonal * self.A_e
        A[:J, 4 * J:5 * J] = displacement_tridiagonal * self.A_v

        # ---
        # Boundary Conditions

        # no displacement at skull
        A[J - 1, J - 1] = 1.0
        b[J - 1] = 0.0

        # constant arterial blood pressure at the skull
        A[2 * J - 1, 2 * J - 1] = 1.0
        b[2 * J - 1] = self.p_bpA

        # no capillary flow at the skull
        A[3 * J - 1, 3 * J - 2] = -1.0
        A[3 * J - 1, 3 * J - 1] = 1.0
        b[3 * J - 1] = 0.0

        # CSF pressure at skull
        A[4 * J - 1, 4 * J - 1] = 1.0
        b[4 * J - 1] = self.p_bp + self.mu_e * self.R * self.Q_o

        # constant venous blood pressure at the skull
        A[5 * J - 1, 5 * J - 1] = 1.0
        b[5 * J - 1] = self.p_bp

        # stress equilibrium in the ventricle wall
        A[0, 0] = (2. * self.nu / r[0] - (1. - self.nu) / dr)
        A[0, 1] = (1. - self.nu) / dr
        b[0] = 0.0
        stress_multiplier = (1. + self.nu) * (1. - 2. * self.nu) / self.E
        A[0, 1 * J] = (1. - self.A_a) * stress_multiplier
        A[0, 2 * J] = (1. - self.A_c) * stress_multiplier
        A[0, 3 * J] = (1. - self.A_e) * stress_multiplier
        A[0, 4 * J] = (1. - self.A_v) * stress_multiplier

        # no arteriol blood flow into ventricles
        A[J, J] = -1.0
        A[J, J + 1] = 1.0
        b[J] = 0.0

        # capillary blood flow into ventricles
        A[2 * J, 2 * J] = -1.0
        A[2 * J, 2 * J + 1] = 1.0
        b[2 * J] = self.mu_a * dr * self.Q_p / self.k_ce

        # venous blood flow into ventricles
        A[4 * J, 4 * J] = -1.0
        A[4 * J, 4 * J + 1] = 1.0
        b[4 * J] = 0.0

        # ---
        # Transient Boundary Conditions
        d2 = self.d ** 2
        d4 = d2 * d2
        const_1 = np.pi * d4 / (128. * self.mu_e * self.L)
        const_2 = 4. * np.pi * (r[0] + x[0]) * (r[0] + x[0])
        # conservation of mass in ventricle
        A[3 * J, 0] = const_2 / self.dt
        A[3 * J, 3 * J] = const_1 + const_2 * self.kappa_e / dr
        A[3 * J, 3 * J + 1] = -const_2 * self.kappa_e / dr
        A[3 * J, 4 * J - 1] = -const_1
        b[3 * J] = self.Q_p + const_2 * x[0] / self.dt

        # ---
        # Scale to better condition the problem
        for i in range(self.b.shape[0]):
            A_max = np.amax(self.A[i, :])
            A_min = np.amin(self.A[i, :])
            scale = 1.0 / (A_max - A_min)
            self.A[i, :] *= scale
            self.b[i] *= scale
        # plt.imshow(self.A, interpolation='nearest')
        # plt.colorbar()
        # plt.show()
        # plt.clf()

    def _save_wall_file(self):
        with open(self.wall_file_name, "w") as f:
            fmt_str = ", ".join(["{:.6e}"] * 11) + os.linesep
            f.write("r, u, p_a, p_c, p_e, p_v, " +
                    "u_res, p_a_res, p_c_res, p_e_res, p_v_res" +
                    os.linesep)
            for i in range(self.J):
                f.write(fmt_str.format(self.r[i],
                                       self.x[i + 0 * self.J],
                                       self.x[i + 1 * self.J],
                                       self.x[i + 2 * self.J],
                                       self.x[i + 3 * self.J],
                                       self.x[i + 4 * self.J],
                                       self.residual[i + 0 * self.J],
                                       self.residual[i + 1 * self.J],
                                       self.residual[i + 2 * self.J],
                                       self.residual[i + 3 * self.J],
                                       self.residual[i + 4 * self.J]))

    def _create_transient_file(self):
        with open(self.transient_fileName, "w") as f:
            f.write("T, U, Pa, Pc, Pe, Pv, " +
                    "U_res, Pa_res, Pc_res, Pe_res, Pv_res" +
                    os.linesep)

    def _save_transient_file(self):
        fmt_str = ", ".join(["{:.6e}"] * 11) + os.linesep
        with open(self.transient_fileName, "a") as f:
            f.write(fmt_str.format(self.t,
                                   self.x[0 * self.J],
                                   self.x[1 * self.J],
                                   self.x[2 * self.J],
                                   self.x[3 * self.J],
                                   self.x[4 * self.J],
                                   self.residual[0 * self.J],
                                   self.residual[1 * self.J],
                                   self.residual[2 * self.J],
                                   self.residual[3 * self.J],
                                   self.residual[4 * self.J]))

    def solve(self):
        if self.save_transient:
            self._create_transient_file()

        for i in range(self.N):
            self.t = self.t0 + i * self.dt

            self._build_system()
            self.x = spsl.spsolve(sps.csc_matrix(self.A), self.b)
            self.residual = np.dot(self.A, self.x) - self.b

            if self.debug_print:
                relative_error = np.linalg.norm(self.residual) / np.linalg.norm(self.b)  # norm() is L2 norm
                print self.base_name,
                print "- Current time:", self.t, " sec. The relative error is:", relative_error

            if self.save_transient:
                self._save_transient_file()

            if self.x[0] > 1e5:
                # simulation is getting too big... cancel solution
                break

        if self.save_wall:
            self._save_wall_file()

