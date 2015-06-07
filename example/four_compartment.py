import numpy as np


class FourCompartmentMPET(object):
    def __init__(self,
                 # general constants
                 final_time, dt, grid_spacing, name,
                 # arteriol constants
                 A_a, B_a, k_a, mu_a,
                 # capillary constants
                 A_c, B_c, k_c, mu_c, k_ce,
                 # venus constants
                 A_v, B_v, k_v, mu_v,
                 # transfer constants
                 gamma_ac, gamma_ce, gamma_cv, gamma_ev,
                 blocked=False):

        # geometry constants
        self.rV = 30.0e-3  # m
        self.rS = 100.0e-3  # m
        self.L = self.rS - self.rV  # m
        self.d = 0.25e-3 if blocked else 4e-3  # m (assume a blocked aqueduct 0.25mm)

        # poroelastic constants
        self.nu = 0.35
        self.E = 584.0  # N/m^2
        self.G = self.E / 2. / (1. + self.nu)  # N/m^2
        self.K = self.E / 3. / (1. - 2. * self.nu)  # N/m^2

        # arteriol constants
        self.A_a = A_a
        self.B_a = B_a
        self.M_a = self.B_a * self.K / self.A_a / (1 - self.B_a * self.A_a)  # N/m^2
        self.k_a = k_a  # m^2
        self.mu_a = mu_a  # Ns/m^2
        self.kappa_a = self.k_a / self.mu_a  # m^4/Ns

        # capillary constants
        self.A_c = A_c
        self.B_c = B_c
        self.M_c = self.B_c * self.K / self.A_c / (1 - self.B_c * self.A_c)  # N/m^2
        self.k_c = k_c  # m^2
        self.mu_c = mu_c  # Ns/m^2
        self.kappa_c = self.k_c / self.mu_c  # m^4/Ns
        self.k_ce = k_ce

        # extracellular/CSF constants
        self.A_e = 1.0
        self.B_e = 0.99
        self.M_e = self.B_e * self.K / self.A_e / (1 - self.B_e * self.A_e)  # N/m^2
        self.k_e = 1.4e-14  # m^2
        self.mu_e = 8.9e-4  # Ns/m^2
        self.kappa_e = self.k_e / self.mu_e  # m^4/Ns

        # venous constants
        self.A_v = A_v
        self.B_v = B_v
        self.M_v = self.B_v * self.K / self.A_v / (1 - self.B_v * self.A_v)  # N/m^2
        self.k_v = k_v  # m^2
        self.mu_v = mu_v  # Ns/m^2
        self.kappa_v = self.k_v / self.mu_v  # m^4/Ns

        # transfer coefficients
        self.gamma_ac = gamma_ac
        self.gamma_ce = gamma_ce
        self.gamma_cv = gamma_cv
        self.gamma_ev = gamma_ev

        # flow contants
        self.Q_p = 5.8e-9  # m^3/s
        self.Q_o = self.Q_p * 1.0  # m^3/s
        self.R = 8.5e13  # m^-3
        self.p_bp = 650.0  # N/m^2
        self.p_bpA = 13.3e3  # N/m^2 arterial blood pressure (100mmHg)

        # grid properties
        self.J = grid_spacing
        self.r = np.linspace(self.rV, self.rS, self.J)
        self.dr = self.r[1] - self.r[0]

        # time control
        self.dt = dt  # time step size
        self.tf = final_time  # final solution time
        self.t0 = 0.0  # initiali solution time
        self.t = self.t0  # current time
        self.N = int((self.tf - self.t0 + self.dt / 10.) / self.dt) + 1  # number of time steps

        # solution properties
        self.x = None  # LHS (solution) vector
        self.b = None  # RHS vector
        self.A = None  # Derivative matrix

        # output file constants
        self.base_name = name
        self.transientFileName = None
        self.wallFileName = None
        self.saveTrans = None

    def _initialise_system(self):

        # == build the A matrix ==
        n = 5  # U, Pa, Pc, Pe, Pv
        vec_size = n * self.J
        self.x = np.zeros((vec_size,))
        self.b = np.zeros((vec_size,))
        self.A = np.zeros((vec_size, vec_size))

        self._initialise_difference_tridiagonals()
        self._initialise_poro_diagonals()
        self._initialise_boundary_conditions()

    def _initialise_difference_tridiagonals(self):
        A = self.A
        J = self.J
        dr = self.dr
        dr2 = dr * dr
        r = self.r
        central_r = self.r[1:-1]
        # set up the difference tridiagonals
        difference_eye_minus = np.diag(1. / dr2 - 1. / (r[:-1] * dr), -1)
        difference_eye_central = np.eye(J, k=0) * (-2.0 / dr2)
        difference_eye_plus = np.diag(1. / dr2 + 1. / (r[1:] * dr), 1)
        difference_tridiagonal = difference_eye_minus + difference_eye_central + difference_eye_plus
        difference_tridiagonal[0, :] = 0.0
        difference_tridiagonal[-1, :] = 0.0
        for i in range(self.x.shape[0] / self.J):
            A[i * J:(i + 1) * J, i * J:(i + 1) * J] = difference_tridiagonal
        # update for the displacement equation
        A[1:J - 1, 1:J - 1] += np.diag(- 2. / (central_r * central_r), 0)

    def _initialise_poro_diagonals(self):
        J = self.J
        dr = self.dr
        A = self.A
        # ---
        # set up the poro diagonals
        stiffness_multiplier = (1.0 - 2.0 * self.nu) / (4.0 * self.G * (1.0 - self.nu) * dr)
        displacement_tridiagonal = (np.eye(J, k=-1) - np.eye(J, k=1)) * stiffness_multiplier
        displacement_tridiagonal[0, :] = 0.0
        displacement_tridiagonal[-1, :] = 0.0
        A[:J, 1 * J:2 * J] = displacement_tridiagonal * self.A_a
        A[:J, 2 * J:3 * J] = displacement_tridiagonal * self.A_c
        A[:J, 3 * J:4 * J] = displacement_tridiagonal * self.A_e
        A[:J, 4 * J:5 * J] = displacement_tridiagonal * self.A_v

    def _initialise_boundary_conditions(self):
        J = self.J
        dr = self.dr
        r = self.r
        A = self.A
        b = self.b

        # ---
        # Boundary Conditions

        # no displacement at skull
        A[J - 1, J - 1] = 1.0

        # constant arterial blood pressure at the skull
        A[2 * J - 1, 2 * J - 1] = 1.0
        b[2 * J - 1] = self.p_bpA

        # no capillary flow at the skull
        A[3 * J - 1, 3 * J - 2] = -1.0
        A[3 * J - 1, 3 * J - 1] = 1.0

        # CSF pressure at skull
        A[4 * J - 1, 4 * J - 1] = 1.0
        b[4 * J - 1] = self.p_bp + self.mu_e * self.R * self.Q_o

        # constant venous blood pressure at the skull
        A[5 * J - 1, 5 * J - 1] = 1.0
        b[5 * J - 1] = self.p_bp

        # stress equilibrium in the ventricle wall
        A[0, 0] = 2. * self.nu / r[0] - (1. - self.nu) / dr
        A[0, 1] = (1. - self.nu) / dr
        stress_multiplier = (1. + self.nu) * (1. - 2. * self.nu) / self.E
        A[0, 1 * J] = (1. - self.A_a) * stress_multiplier
        A[0, 2 * J] = (1. - self.A_c) * stress_multiplier
        A[0, 3 * J] = (1. - self.A_e) * stress_multiplier
        A[0, 4 * J] = (1. - self.A_v) * stress_multiplier

        # no arteriol blood flow into ventricles
        A[J, J] = -1.0
        A[J, J + 1] = 1.0

        # capillary blood flow into ventricles
        A[2 * J, 2 * J] = -1.0
        A[2 * J, 2 * J + 1] = 1.0
        b[2 * J] = self.mu_a * dr * self.Q_p / self.k_ce

        # venous blood flow into ventricles
        A[4 * J, 4 * J] = -1.0
        A[4 * J, 4 * J + 1] = 1.0

    def _update_time_dependent_system(self):
        # == build the A matrix ==
        J = self.J
        dr = self.dr
        r = self.r
        x = self.x
        b = self.b
        A = self.A

        # ---
        # ** Pressure equations
        sDot_ac = self.gamma_ac * np.abs(x[0 * J + 1:1 * J - 1] - x[1 * J + 1:2 * J - 1])
        sDot_ce = self.gamma_ce * np.abs(x[2 * J + 1:3 * J - 1] - x[3 * J + 1:4 * J - 1])
        sDot_cv = self.gamma_cv * np.abs(x[2 * J + 1:3 * J - 1] - x[4 * J + 1:5 * J - 1])
        sDot_ev = self.gamma_ev * np.abs(x[3 * J + 1:4 * J - 1] - x[4 * J + 1:5 * J - 1])
        # arteriol pressure equation
        b[J + 1:2 * J - 1] = sDot_ac / self.kappa_a
        # capillary pressure equation
        b[2 * J + 1:3 * J - 1] = (-sDot_ac + sDot_ce + sDot_cv) / self.kappa_c
        # CSF pressure equation
        b[3 * J + 1:4 * J - 1] = (-sDot_ce + sDot_ev) / self.kappa_e
        # venous pressure equation
        b[4 * J + 1:5 * J - 1] = (-sDot_cv - sDot_ev) / self.kappa_v

        # ---
        # Boundary Conditions
        d4 = self.d ** 2 * self.d ** 2
        const_1 = np.pi * d4 / (128. * self.mu_e * self.L)
        const_2 = 4. * np.pi * np.power(r[0] + x[0], 2)
        # conservation of mass in ventricle
        A[3 * J, 0] = const_2 / self.dt
        A[3 * J, 3 * J] = const_1 + const_2 * self.kappa_e / dr
        A[3 * J, 3 * J + 1] = -const_2 * self.kappa_e / dr
        A[3 * J, 4 * J - 1] = -const_1
        b[3 * J] = self.Q_p + const_2 * x[0] / self.dt

    def solve(self):
        self._initialise_system()
        for i in range(self.N):
            self.t = self.t0 + i * self.dt
            self._update_time_dependent_system()
            self.x = np.linalg.solve(self.A, self.b)
            if self.x[0] > 1e5:
                # simulation is getting too big... cancel solution
                break


def main():
    grid_spacing = 81
    secs_in_day = 86400
    final_time = 1.0 * secs_in_day
    dt = 100.0
    name = "example"

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

    solver = FourCompartmentMPET(final_time, dt, grid_spacing, name,
                                 alpha_a, beta_a, k_a, mu_a,
                                 alpha_c, beta_c, k_c, mu_c, k_ce,
                                 alpha_v, beta_v, k_v, mu_v,
                                 gamma_ac, gamma_ce, gamma_cv, gamma_ev)
    solver.solve()
    # import matplotlib.pyplot as plt
    # plt.plot(solver.r, solver.x[:solver.J])
    # plt.show()
    # print solver.x


if __name__ == "__main__":
    main()
