from mpet import FourCompartmentPoro


def main():
    grid_spacing = 151
    secs_in_day = 86400
    initial_time = 0.0
    final_time = 1.0 * secs_in_day
    dt = 100.0
    write_transient = False
    base_name = "example"

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

    s = FourCompartmentPoro()

    s.initialize(grid_spacing, initial_time, final_time, dt, write_transient, base_name)
    s.setArteriolConstants(alpha_a, beta_a, k_a, mu_a)
    s.setCapillaryConstants(alpha_c, beta_c, k_c, mu_c, k_ce)
    s.setVenousConstants(alpha_v, beta_v, k_v, mu_v)
    s.setTransferConstants(gamma_ac, gamma_ce, gamma_cv, gamma_ev)
    s.solve()


if __name__ == "__main__":
    main()
