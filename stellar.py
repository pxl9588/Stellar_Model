from math import *
try:
    from tabulate import *
except:
    print("pip install tabulate brah")
GAMMA = 5/3
constants = {"sigma": 5.67051e-5, "c": 2.99792458e+10, "a": 7.56591e-15, "G": 6.67259e-8, "k_B": 1.38065e-16,
                 "m_H": 1.673534e-24, "pi": 3.141592654, "gamma":  GAMMA, "gamrat":  GAMMA / (GAMMA - 1.0),
                 "kPad":  1.0, "g_ff": 1.0, "Rsun": 6.9599e+10, "Msun": 1.989e+33, "Lsun": 3.826e+33}
diff_eqs = {"dMdr": 0, "dPdr": 0, "dLdr": 0, "dTdr": 0}

def statstar():
    #  FLAGS:
    #  ---------------------------------------------
    #  idrflg {0 = Rs/1000, 1 = Rs/100 , 2 = Rs/5000} initial dr flag
    #  i_goof {-1,0,1,2,3,4,5}
    #  ----------------------------------------------
    flags = {"idrflg": 0, "i_goof": -1}

    # Goof Dictionary
    goof = {1: "1"}

    # N_start number of steps for which starting equations are used
    # N_stop maximum number of allowed zones
    zone_boundaries = {"N_start": 10, "N_stop": 999}  # N_start,N_stop

    #kPad = p/t^gamrat

    # initial_cond = [Msolar, Lsolar, Te]
    initial_cond = {"Msolar": 0, "Lsolar": 0, "Rsolar": 0, "Te": 0, "Ms": 0, "Ls": 0, "Rs": 0, "T0": 0, "P0": 0}

    # mass_fractions = [X, Y, Z]
    mass_fractions = {"X": 0, "Y": 0, "Z": 0, "XCNO": 0}

    # OPEN SOME FILE FOR SOMETHING

    initial_cond["Msolar"] = 1.0  # float(input("Enter the mass of the star (in solar units): "))

    initial_cond["Lsolar"] = 1.0  # float(input("Enter the luminosity of the star (in solar units): "))

    initial_cond["Te"] = 5600.0  # float(input("Enter the effective temperature of the star (in K): "))

    mass_fractions["X"] = .75  # float(input("Enter the mass fraction of hydrogen (X): "))

    mass_fractions["Z"] = .02  # float(input("Enter the mass fraction of metals (Z): "))

    mass_fractions["Y"] = 1.0 - mass_fractions["X"] - mass_fractions["Z"]
    if mass_fractions["Y"] < 0.0:
        SystemExit("Where's the hydrogen?")

    # Select the mass fraction CNO to be 50% of Z.
    mass_fractions["XCNO"] = mass_fractions["Z"] / 2.0

    # Calculate the mass, luminosity, and radius of the star
    initial_cond["Ms"] = initial_cond["Msolar"] * constants["Msun"]
    initial_cond["Ls"] = initial_cond["Lsolar"] * constants["Lsun"]
    initial_cond["Rs"] = ((initial_cond["Ls"]/(4.0*constants["pi"]*constants["sigma"]))**1/2)/(initial_cond["Te"]**2)
    initial_cond["Rsolar"] = initial_cond["Rs"] / constants["Rsun"]

    resultDictionary = {"zone_boundaries": zone_boundaries, "initial_cond": initial_cond,
                        "mass_fractions": mass_fractions, "constants": constants}

    # Begin with a very small step size since surface conditions vary rapidly.
    delta_r = -initial_cond["Rs"] / 1000.0

    # Calculate mean molecular weight mu assuming complete ionization
    mu = 1.0 / (2.0 * mass_fractions["X"] + .75 * mass_fractions["Y"] + .5 * mass_fractions["Z"])

    # Initialize values of r, p, m_r, l_r, t, rho, kappa, and epsilon at the surface.
    r, m_r, l_r, t, p, rho, kappa, epsilon, dl_pdl_t = ([0.0] * zone_boundaries["N_stop"] for makeArray in range(9))

    # took tog_bof out of constants because its dynamic
    tog_bf = .01
    dl_p_lim = 99.9

    r[0] = initial_cond["Rs"]
    m_r[0] = initial_cond["Ms"]
    l_r[0] = initial_cond["Ls"]
    t[0] = initial_cond["T0"]
    p[0] = initial_cond["P0"]
    if p[0] <= 0.0 or t[0] <= 0.0:
        rho[0] = 0.0
        kappa[0] = 0.0
        epsilon[0] = 0.0
    else:
        [rho[0], kappa[0], epsilon[0], tog_bf] = eos(mass_fractions["X"], mass_fractions["Z"], mass_fractions["XCNO"], mu, p[0], t[0])

    constants["kPad"] = 0.3
    irc = 0
    dl_pdl_t[0] = 4.25
    rho_core, eps_core, p_core, t_core, i_stop, i_goof = 0, 0, 0, 0, 0, 0

    # do 20 i = 1, Nstart
    for i in range(0, zone_boundaries["N_start"] + 1):
        ip1 = i + 1
        [r[ip1], m_r[ip1], l_r[ip1], t[ip1], p[ip1]] = startmdl(delta_r, mass_fractions["X"], mass_fractions["Z"], mu, initial_cond["Rs"], r[i], m_r[i], l_r[i], tog_bf, irc)
        [rho[ip1], kappa[ip1], epsilon[ip1], tog_bf] = eos(mass_fractions["X"], mass_fractions["Z"], mass_fractions["XCNO"], mu, p[ip1], t[ip1])

        # Determine whether convection will be operating in th next zone
        if i > 1:
            dl_pdl_t[ip1] = log(p[ip1]/p[i])/log(t[ip1]/t[i])
        else:
            dl_pdl_t[ip1] = dl_pdl_t[i]
        if dl_pdl_t[ip1] < constants["gamrat"]:
            irc = 1
        else:
            irc = 0
            constants["kPad"] = p[ip1]/t[ip1]**constants["gamrat"]

        # Test to see whether the surface assumption of constant mass is still valid
        delta_m = delta_r * dm_dr(r[ip1], rho[ip1])
        m_r[ip1] = m_r[i] + delta_m
        if abs(delta_m) > .001 * initial_cond["Ms"]:
            print("The variation in mass has become larger than 0.001*Ms, leaving the approximation loop before Nstart was reached")
            if ip1 > 2:
                ip1 = ip1 - 1
            break

    # Main integration loop
    nstrtip1 = ip1 + 1
    f_im1 = [0] * 4
    df_dr = [0] * 4
    f_i = [0, 0, 0, 0]

    for j in range(nstrtip1, zone_boundaries["N_stop"]):
        im1 = j - 1

        # Initialize the Runge-Kutta routine with zone i - 1 quantities and their derivatives.
        # The pressure, mass, luminosity and temperature are stored in f_im(0-3).
        # The derivatives of those quantities with respect to radius are in df_dr(0-3).
        # The The resulting values for p, m_r, l_r, and t are returned in f_i(0-3)

        f_im1[0] = p[im1]
        f_im1[1] = m_r[im1]
        f_im1[2] = l_r[im1]
        f_im1[3] = t[im1]
        df_dr[0] = dp_dr(r[im1], m_r[im1], rho[im1])
        df_dr[1] = dm_dr(r[im1], rho[im1])
        df_dr[2] = dl_dr(r[im1], rho[im1], epsilon[im1])
        df_dr[3] = dt_dr(r[im1], m_r[im1], l_r[im1], t[im1], rho[im1], kappa[im1], mu, irc)
        runge(f_im1, df_dr, f_i, r[im1], delta_r, irc, mass_fractions["X"], mass_fractions["Z"], mass_fractions["XCNO"], mu)

        # Update stellar parameters for the next zone, including adding dr to the old radius

        r[j] = r[im1] + delta_r
        p[j] = f_i[0]
        m_r[j] = f_i[1]
        l_r[j] = f_i[2]
        t[j] = f_i[3]

        # Calculate the density, opacity, and energy generation rate for this zone
        [rho[j], kappa[j], epsilon[j], tog_bf, ierr] = eos(mass_fractions["X"], mass_fractions["Z"], mass_fractions["XCNO"], mu, p[j], t[j])

        # Determine whether convection will be operating in the next zone
        dl_pdl_t[j] = log(p[j]/p[im1])/log(t[j]/t[im1])
        if dl_pdl_t[j] < constants["gamrat"]:
            irc = 1
        else:
            irc = 0

        # Check if the center has been reached. If so, set i_goof and estimate the central conditions:
        # rho_core, eps_core, p_core, t_core.

        if r[j] < abs(delta_r) and (l_r[j] > .1 * initial_cond["Ls"] or m_r[j] > .01 * initial_cond["Ms"]):
            i_goof = 6
        elif l_r[j] < 0.0:
            i_goof = 5
            rho_core = m_r[j]/(4/3*pi*r[j]**3)
            if m_r[j] != 0:
                eps_core = l_r[j]/m_r[j]
            else:
                eps_core = 0.0
            p_core = p[j] + 2/3 * pi * constants["G"] * rho_core**2 * r[j]**2
            t_core = p_core * mu * constants["m_H"]/(rho_core*constants["k_B"])
        elif m_r[j] < 0.0:
            i_goof = 4
            rho_core = 0.0
            eps_core = 0.0
            p_core = 0.0
            t_core = 0.0
        elif r[j] < .02 * initial_cond["Rs"] and m_r[j] < .01 * initial_cond["Ms"] and l_r[j] < .1 * initial_cond["Ls"]:
            rho_core = m_r[j] / (4/3 * pi * r[j]**3)
            rho_max = 10.0 * (rho[j]/rho[im1]) * rho[j]
            eps_core = l_r[j] / m_r[j]
            p_core = p[j] + 2/3 * pi * constants["G"] * rho_core**2 * r[j]**2
            t_core = p_core * mu * constants["m_H"] / (rho_core * constants["k_B"])
            if rho_core < rho[j] or rho_core > rho_max:
                i_goof = 1
            elif eps_core < epsilon[j]:
                i_goof = 2
            elif t_core < t[j]:
                i_goof = 3
            else:
                i_goof = 0
        if i_goof != -1:
            i_stop = j
            # Loop
            break

        # Change step size?
        if flags["idrflg"] == 0 and m_r[j] < .99 * initial_cond["Ms"]:
            delta_r = -initial_cond["Rs"]/100.0
            flags["idrflg"] = 1
        if flags["idrflg"] == 1 and delta_r > .5 * r[j]:
            delta_r = -initial_cond["Rs"]/5000
            flags["idrflg"] = 2
        i_stop = j
    # End of loop

    rho_core = m_r[i_stop] / (4/3 * pi * r[i_stop]**3)
    eps_core = l_r[i_stop] / m_r[i_stop]
    p_core = p[i_stop] + 2 / 3 * pi * constants["G"] * rho_core ** 2 * r[i_stop] ** 2
    t_core = p_core * mu * constants["m_H"] / (rho_core * constants["k_B"])

    if i_goof != 0:
        print(goof[i_goof])
    else:
        print("Success")

    # Print the central conditions
    r_crat = r[i_stop] / initial_cond["Rs"]
    if r_crat < -9.999:
        r_crat = -9.999
    m_crat = m_r[i_stop] / initial_cond["Ms"]
    if m_crat < -9.999:
        m_crat = -9.999
    l_crat = l_r[i_stop] / initial_cond["Ls"]
    if l_crat < -9.999:
        l_crat = -9.999

    solar_crat_table = [["", "Solar", "Crat"], ["Mass", initial_cond["Msolar"], m_crat], ["Radius", initial_cond["Rsolar"], r_crat], ["Luminosity", initial_cond["Lsolar"], l_crat]]
    print(tabulate(solar_crat_table), headers="firstrow")
    mass_table = [["Hydrogen", "Helium", "Metals"], ["Mass Fractions", mass_fractions["X"], mass_fractions["Y"], mass_fractions["Z"]]]
    print(tabulate(mass_table), headers="firstrow")
    core_table = [["Rho", "Temperature", "Pressure", "Epsilon"], [rho_core, t_core, p_core, eps_core]]
    print(tabulate(core_table, headers="firstrow"))

    # Print data from the cetner of the start outward, labeling convective or radiative zones zones.
    for ic in range(0, i_stop):
        i = i_stop - ic + 1
        qm = 1.0 - m_r[i] / initial_cond["Ms"]
        '''if dl_pdl_t[i] < constants["gamrat"]:
            rcf = 'c'
        else:
            rcf = 'r'
        if abs(dl_pdl_t[i]) > dl_p_lim:
            if dl_pdl_t[i] >= 0:
                dl_pdl_t = dl_p_lim
            else:
                dl_pdl_t = -dl_p_lim
            clim = '*'
        else:
            clim = ' '
        '''
        results = [["R [" + str(i) + "]", r[i]], ["Qm", qm], ["l_r[" + str(i) + "]", l_r[i]], ["t[" + str(i) + "", t[i]],
                   ["p["+ str(i) + "]", p[i]], ["Rho[" + str(i) + "]", rho[i]], ["Kappa[" + str(i) + "]", kappa[i]],
                   ["Epsilon[" + str(i) + "]", epsilon[i]], ["dl_pdl_t[" + str(i) + "]", dl_pdl_t[i]]]
        print(tabulate(results))
    print("**** The integration has been completed ****")

#calculates the values of density,opacity, guillotine-to-gaunt ration, energy gen rate
# for a radius r

# 'Hydrostatic equilibrium Pressure Gradient
def dp_dr(r, m_r, rho):
    #dPdr
    return -constants["G"] * rho * m_r * r**-2

# 'Conservation of Mass
def dm_dr(r, rho):
    #dMdr
    return 4 * pi * rho * r**2


# 'luminosity thingy
def dl_dr(r, rho, epslon):
    #dLdr
    return 4 * pi * rho * epslon * r**2


# 'The temp one
def dt_dr(r, m_r, l_r, t, rho, kappa, mu, irc):
    if irc == 0:
        #dTdr
        return -(3/(16 * pi * constants["a"] * constants["c"])) * kappa * rho * t ** -3 * l_r * r ** -2
    else:
        #dTdr
        return -1 / constants["gamrat"] * constants["G"] * m_r * r ** -2 * mu * constants["m_H"] / constants["k_B"]


# Returns r, M_rip1, L_rip1, T_ip1, P_ip1
def startmdl(delta_r, x, z, mu, rs, r_i, m_ri, l_ri, tog_bf, irc):

    r = r_i + delta_r
    m_rip1 = m_ri
    l_rip1 = l_ri

    if irc == 0:
        t_ip1 = constants["G"] * m_rip1 * mu * constants["m_H"]/(4.25*constants["k_B"])*(1.0 / r - 1.0 / rs)
        a_bf = 4.34e+25 * z * (1.0 + x) / tog_bf
        a_ff = 3.68e+22 * constants["g_ff"] * (1.0 - z) * (1.0 + x)
        afac = a_bf + a_ff

        p_ip1 = sqrt((1.0/4.25)*(16.0/3.0*pi*constants["a"]*constants["c"])*(constants["G"] * m_rip1/l_rip1) *
                     (constants["k_B"]/(afac * mu * constants["m_H"]))) * t_ip1**4.25

    # This is the convective approximation
    else:
        t_ip1 = constants["G"] * m_rip1 * mu * constants["m_H"] / constants["k_B"] * (1.0 / r - 1.0 / rs) / constants["gamrat"]
        p_ip1 = constants["kPad"]*t_ip1**constants["gamrat"]

    return [r, m_rip1, l_rip1, t_ip1, p_ip1]

def runge(f_im1, dfdr, f_i, r_im1, deltar, irc, x, z, xcno, mu):
    f_temp = [0] * 4
    df1 = [0] * 4
    df2 = [0] * 4
    df3 = [0] * 4
    dr12 = deltar/2.0
    dr16 = deltar/6.0
    r12 = r_im1 + dr12
    r_i = r_im1 + deltar
    # Calculate intermediate derivatives from the fundamental stellar structure equations found in subroutine fundeq.
    for i in range(0, 4):
        f_temp[i] = f_im1[i] + dr12*dfdr[i]
    fundeq(r12, f_temp, df1, irc, x, z, xcno, mu)
    for i in range(0, 4):
        f_temp[i] = f_im1[i] + dr12*df1[i]
    fundeq(r12, f_temp, df2, irc, x, z, xcno, mu)
    for i in range(0, 4):
        f_temp[i] = f_im1[i] + deltar*df2[i]
    fundeq(r_i, f_temp, df3, irc, x, z, xcno, mu)
    for i in range(0, 4):
        f_i[i] = f_im1[i] + dr16*(dfdr[i] + 2*df1[i] + 2.0 * df2[i] + df3[i])


def fundeq(r, f, dfdr, irc, x, z, xcno, mu):
    p = f[0]
    m_r = f[1]
    l_r = f[2]
    t = f[3]
    [rho, kappa, epslon, tog_bf] = eos(x, z, xcno, mu, p, t)
    dfdr[0] = dp_dr(r, m_r, rho)
    dfdr[1] = dm_dr(r, rho)
    dfdr[2] = dl_dr(r, rho, epslon)
    dfdr[3] = dt_dr(r, m_r, l_r, t, rho, kappa, mu, irc)


# instead of passing in ierr going to use the return as an error code. 0 for success.
# Returns [rho, kappa, epslon, tog_bf]
def eos(x, z, xcno, mu, p, t):

    oneo3 = .333333333333
    twoo3 = .666666666667

    # solve for density from the ideal gas law
    if t <= 0.0 or p <= 0.0:
        raise SystemExit('Temperature or Pressure less than equal to 0, in EOS')
    prad = constants['a'] * (t ** 4) / 3.0
    pgas = p - prad
    rho = (mu * constants['m_H']/constants['k_B']) * (pgas / t)
    if rho < 0.0:
        raise SystemExit("Rho is less than 0, in EOS")
    # Calc opacity, including guillatine-to-gaunt factor ratio

    tog_bf = 2.82 * (rho * (1 + x)) ** .2
    k_bf = 4.34E25 / tog_bf * z(1 + x) * rho / (t ** 3.5)
    k_ff = 3.68E22 * constants["g_ff"] * (1 - z) * (1 + x) * rho / (t ** 3.5)
    k_e = .2*(1 + x)
    kappa = k_bf + k_ff + k_e

    #calc eneregy generation by pp chain and CNO cycle
    #The screening factor for the pp chain is calculated as fpp

    t6 = t * 1E-6
    fx = .133 * x * ((3 + x) * rho) ** .5 / t6 ** 1.5
    fpp = 1 + fx * x
    psipp = 1 + 1.412E8 * (1 / x - 1) * exp(-49.98 * t6 ** -oneo3)
    cpp = 1 + .0123*t6*oneo3 + .0109*t6**(-twoo3) - .000149*t6
    epspp = 2.38E6 * rho * x * x * fpp * psipp * cpp * t6 ** (-twoo3) * exp(-33.8 * t6 ** (-oneo3))
    ccno = 1 + .0027*t6**oneo3 - .0078 * t6**twoo3 - .000149*t6
    eps_cno = 8.67E27 * rho * x * xcno * ccno * t6 ** (-twoo3) * exp(-152.28 * t6 ** (-oneo3))
    epslon = epspp + eps_cno

    return [rho, kappa, epslon, tog_bf]


statstar()