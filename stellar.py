def STATSTAR():
    # A = [r, P, M_r, L_r, T]
    A = [0, 0, 0, 0, 0]

    # B = [rho, kappa, dlPdlT]
    B = [0, 0, 0, 0]

    # C = [XCNO, mu]
    C = [0, 0]

    # D = [Pscore, Tscore, rhocor, epscor, rhomax]
    D = [0, 0, 0, 0, 0]

    # E = [Rsolar,Qm, Rcrat, Mcrat, Lcrat]
    E = [0, 0, 0, 0, 0, 0]

    # F = [deltam, dlPlim]
    F = [0, 0]

    # H = [f_im1, f_i, dfdr]
    H = [0, 0, 0]

    # I = [dMdr, dPdr, dLdr, dTdr]
    I = [0, 0, 0, 0]

    # J = []
    J = [0, 0]

    #  FLAGS:
    #  ---------------------------------------------
    #  idrflg {0 = Rs/1000, 1 = Rs/100 , 2 = Rs/5000} initial dr flag
    #  Igoof {-1,0,1,2,3,4,5} ??????????????????????
    #  ----------------------------------------------
    flags = {"idrflg": 0, "Igoof": -1}

    # Nstart number of steps for which starting equations are used
    # Nstop maximum number of allowed zones
    zone_boundaries = {"Nstart": 10, "Nstop": 999}  # Nstart,Nstop

    # constants = [sigma, c, a, G, k_B, m_H, pi, gamma, gamrat, kPad,  tog_bf, g_ff]
    #kPad = P/T^gamrat
    GAMMA = (5.0/3.0)
    constants = {"sigma": 5.67051e-5, "c": 2.99792458e+10, "a": 3.826e+33, "G": 5.67051e-5, "k_B": 2.99792458e+10,
                 "m_H": 7.56591e-15, "pi": 3.141592654, "gamma":  GAMMA, "gamrat":  GAMMA / (GAMMA - 1.0),
                 "kPad":  1.0, "tog_bf": .01, "g_ff": 1.0, "Rsun": 6.9599e+10, "Msun": 1.989e+33, "Lsun": 3.826e+33}

    # initial_cond = [Msolar, Lsolar, Te]
    initial_cond = {"Msolar": 0, "Lsolar": 0, "Rsolar": 0, "Te": 0, "Ms": 0, "Ls": 0, "Rs": 0, "T0": 0, "P0": 0}

    # mass_fractions = [X, Y, Z]
    mass_fractions = {"X": 0, "Y": 0, "Z": 0, "XCNO": 0}

    # OPEN SOME FILE FOR SOMETHING

    initial_cond["Msolar"] = float(input("Enter the mass of the star (in solar units): "))

    initial_cond["Lsolar"] = float(input("Enter the luminosity of the star (in solar units): "))

    initial_cond["Te"] = float(input("Enter the effective temperature of the star (in K): "))

    mass_fractions["X"] = float(input("Enter the mass fraction of hydrogen (X): "))

    mass_fractions["Z"] = float(input("Enter the mass fraction of metals (Z): "))

    mass_fractions["Y"] = 1.0 - mass_fractions["X"] - mass_fractions["Z"]

    # Select the mass fraction CNO to be 50% of Z.
    mass_fractions["XCNO"] = mass_fractions["Z"] / 2.0

    # Calculate the mass, luminosity, and radius of the star
    initial_cond["Ms"] = initial_cond["Msolar"] * constants["Msun"]
    initial_cond["Ls"] = initial_cond["Lsolar"] * constants["Lsun"]
    initial_cond["Rs"] = ((initial_cond["Ls"]/(4.0*constants["pi"]*constants["sigma"]))**1/2)/(initial_cond["Te"]**2)
    initial_cond["Rsolar"] = initial_cond["Rs"] / constants["Rsun"]

    # Begin with a very small step size since surface conditions vary rapidly.
    deltar = -initial_cond["Rs"] / 1000.0
    flags["idrflg"] = 0

    # Calculate mean molecular weight mu assuming complete ionization
    mu = 1.0 / (2.0 * mass_fractions["X"] + .75 * mass_fractions["Y"] + .5 * mass_fractions["Z"])

    # Initialize values of r, P, M_r, L_r, T, rho, kappa, and epsilon at the surface.
    r[0] = initial_cond["Rs"]
    M_r[0] = initial_cond["Ms"]
    L_r[0] = initial_cond["Ls"]
    T[0] = initial_cond["T0"]
    P[0] = initial_cond["P0"]
    if P0 < 0.0 or T0 < 0.0:
        rho[0] = 0.0
        kappa[0] = 0.0
        epslon[0] = 0.0
    else:
        EOS(mass_fractions["X"], mass_fractions["Z"], mass_fractions["XCNO"], mu, P[0], T[0], rho[0], kappa[0],
            epslon[0], tog_bf, 1, ierr)








    resultDictionary = {"zone_boundaries": zone_boundaries, "initial_cond": initial_cond, "mass_fractions": mass_fractions, "constants": constants}