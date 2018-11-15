from math import *
try:
    from tabulate import *
except:
    print("pip install tabulate brah")
GAMMA = 5/3
constants = {"sigma": 5.67051e-5, "c": 2.99792458e+10, "a": 3.826e+33, "G": 5.67051e-5, "k_B": 2.99792458e+10,
                 "m_H": 7.56591e-15, "pi": 3.141592654, "gamma":  GAMMA, "gamrat":  GAMMA / (GAMMA - 1.0),
                 "kPad":  1.0, "tog_bf": .01, "g_ff": 1.0, "Rsun": 6.9599e+10, "Msun": 1.989e+33, "Lsun": 3.826e+33}
diff_eqs = {"dMdr": 0, "dPdr": 0, "dLdr": 0, "dTdr": 0}

def STATSTAR():
    ierr = 1
    # A = [r, P, M_r, L_r, T]
    A = [0, 0, 0, 0, 0]

    # B = [rho, kappa]
    B = [0, 0]

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


    #  FLAGS:
    #  ---------------------------------------------
    #  idrflg {0 = Rs/1000, 1 = Rs/100 , 2 = Rs/5000} initial dr flag
    #  Igoof {-1,0,1,2,3,4,5}
    #  ----------------------------------------------
    flags = {"idrflg": 0, "Igoof": -1}

    # Goof Dictionary
    goof = {1: "1"}

    # Nstart number of steps for which starting equations are used
    # Nstop maximum number of allowed zones
    zone_boundaries = {"Nstart": 10, "Nstop": 999}  # Nstart,Nstop

    #kPad = P/T^gamrat

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

    resultDictionary = {"zone_boundaries": zone_boundaries, "initial_cond": initial_cond,
                        "mass_fractions": mass_fractions, "constants": constants}

    # Begin with a very small step size since surface conditions vary rapidly.
    deltar = -initial_cond["Rs"] / 1000.0
    flags["idrflg"] = 0

    # Calculate mean molecular weight mu assuming complete ionization
    mu = 1.0 / (2.0 * mass_fractions["X"] + .75 * mass_fractions["Y"] + .5 * mass_fractions["Z"])

    # Initialize values of r, P, M_r, L_r, T, rho, kappa, and epsilon at the surface.
    r = [1]
    M_r = [1]
    L_r = [1]
    T = [1]
    P = [1]
    rho = [1]
    kappa = [1]
    epslon = [1]
    dlPdlT = [1]

    r[0] = initial_cond["Rs"]
    M_r[0] = initial_cond["Ms"]
    L_r[0] = initial_cond["Ls"]
    T[0] = initial_cond["T0"]
    P[0] = initial_cond["P0"]
    if initial_cond["P0"] < 0.0 or initial_cond["T0"] < 0.0:
        rho[0] = 0.0
        kappa[0] = 0.0
        epslon[0] = 0.0
    else:
        EOS(mass_fractions["X"], mass_fractions["Z"], mass_fractions["XCNO"], mu, P[0], T[0], rho[0], kappa[0],
            epslon[0], constants["tog_bf"], 1, ierr)

    constants["kPad"] = 0.3
    irc = 0
    dlPdlT[0] = 4.25
    ip1 = 0

    for i in range(0,zone_boundaries["Nstart"]):
        ip1 = i + 1
        STARTMDL(deltar, mass_fractions["X"], mass_fractions["Z"], mu, initial_cond["Rs"], r[i], M_r[i], L_r[i],
                 [r[ip1], M_r[ip1], L_r[ip1], T[ip1], P[ip1]], constants["tog_bf"], irc)
        EOS(mass_fractions["X"], mass_fractions["Z"], mass_fractions["XCNO"], mu, P[ip1], T[ip1], rho[ip1], kappa[ip1],
            epslon[ip1], constants["tog_bf"], ip1, ierr)

        # Determine whether convection will be operating in th next zone
        if i > 1:
            dlPdlT[ip1] = log(P[ip1]/P[i])/log(T[ip1]/T[i])
        else:
            dlPdlT[ip1] = dlPdlT[i]
        if dlPdlT[ip1] < constants["gamrat"]:
            irc = 1
        else:
            irc = 0
            constants["kPad"] = P[ip1]/T[ip1]**constants["gamrat"]

        # Test to see whether the surface assumption of constant mass is still valid
        deltaM = deltar * dMdr(r[ip1], rho[ip1])
        M_r[ip1] = M_r[i] + deltaM
        if abs(deltaM) > .001 * initial_cond["Ms"]:
            if ip1 > 2:
                ip1 = ip1 - 1
            else:
                continue

        # Main integration loop
        Nsrtp1 = ip1 + 1
        f_im1 = []
        dfdr = []
        f_i = []


        for i in range(Nsrtp1, zone_boundaries["Nstop"]):
            im1 = i - 1

            # Initialize the Runge-Kutta routine with zone i - 1 quantites and their derivatives.
            # The pressure, mass, luminosity and temperature are stored in f_im(0-3).
            # The derivatives of those quantites with respect to radis are in dfdr(0-3).
            # The The resulting values for P, M_r, L_r, and T are returned in f_i(0-3)

            f_im1[0] = P[im1]
            f_im1[1] = M_r[im1]
            f_im1[2] = L_r[im1]
            f_im1[3] = T[im1]
            dfdr[0] = dPdr(r[im1], M_r[im1], rho[im1])
            dfdr[1] = dMdr(r[im1], rho[im1])
            dfdr[2] = dLdr(r[im1], rho[im1], epslon[im1])
            dfdr[3] = dTdr(r[im1], M_r[im1], L_r[im1], T[im1], rho[im1], kappa[im1], mu, irc)
            RUNGE(f_im1, dfdr, f_i, r[im1], deltar, irc, mass_fractions["X"], mass_fractions["Z"], mass_fractions["XCNO"], mu, i, ierr)

            # Update stellar parameters for the next zone, including adding dr to the old radius

            r[i] = r[im1] + deltar
            P[i] = f_i[0]
            M_r[i] = f_i[1]
            L_r[i] = f_i[2]
            T[i] = f_i[3]

            # Calculate the density, opacity, and energy generation rate for this zone
            EOS(mass_fractions["X"], mass_fractions["Z"], mass_fractions["XCNO"], mu, P[i], T[i], rho[i], kappa[i],
                epslon[i], constants["tog_bf"], i, ierr)

            # Determine whether convection will be operating in the next zone

            dlPdlT[i] = log(P[i]/P[im1])/log(T[i]/T[im1])
            if dlPdlT[i] < constants["gamrat"]:
                irc = 1
            else:
                irc = 0

            # Check if the center has been reached. If so, set Igoof and estimate the central conditions:
            # rhocor, epscor, Pcore, Tcore.

            if r[i] < abs(deltar) and (L_r[i] > .1 * initial_cond["Ls"] or M_r[i] > .01 * initial_cond["Ms"]):
                Igoof = 6
            elif L_r[i] < 0.0:
                Igoof = 5
                rhocor = M_r[i]/(4/3*pi*r[i]**3)
                if M_r[i] != 0:
                    epscor = L_r[i]/M_r[i]
                else:
                    epscor = 0.0
                Pcore = P[i] + 2/3 * pi * constants["G"] * rhocor**2 * r[i]**2
                Tcore = Pcore * mu * constants["m_H"]/(rhocor*constants["k_B"])
            elif M_r[i] < 0.0:
                Igoof = 4
                rhocor = 0.0
                epscor = 0.0
                Pcore = 0.0
                Tcore = 0.0
            elif r[i] < .02 * initial_cond["Rs"] and M_r[i] < .01 * initial_cond["Ms"] and L_r[i] < .1 * initial_cond["Ls"]:
                rhocor = M_r[i] / (4/3 * pi * r[i]**3)
                rhomax = 10.0 * (rho[i]/rho[im1]) * rho[i]
                epscor = L_r[i] / M_r[i]
                Pcore = P[i] + 2/3 * pi * constants["G"] * rhocor**2 * r[i]**2
                Tcore = Pcore * mu * constants["m_H"] / (rhocor * constants["k_B"])
                if rhocor < rho[i] or rhocor > rhomax:
                    Igoof = 1
                elif epscor < epslon[i]:
                    Igoof = 2
                elif Tcore < T[i]:
                    Igoof = 3
                else:
                    Igoof = 0
            if Igoof != -1:
                istop = i
                # Loop
                continue

            # Change step size?
            if flags["idrflg"] == 0 and M_r[i] < .99 * initial_cond["Ms"]:
                deltar = -initial_cond["Rs"]/100.0
                flags["idrflg"] = 1
            if flags["idrflg"] == 1 and deltar > .5 * r[i]:
                deltar = -initial_cond["Rs"]/5000
                flags["idrflg"] = 2
            istop = i
        # End of loop

        rhocor = M_r[istop] / (4/3 * pi * r[istop]**3)
        epscor = L_r[istop] / M_r[istop]
        Pcore = P[istop] + 2 / 3 * pi * constants["G"] * rhocor ** 2 * r[istop] ** 2
        Tcore = Pcore * mu * constants["m_H"] / (rhocor * constants["k_B"])

        if Igoof != 0:
            print(goof[Igoof])
        else:
            print("Success")

    # Print the central conditions
    Rcrat = r[istop] / initial_cond["Rs"]
    if Rcrat < -9.999:
        Rcrat = -9.999
    Mcrat = M_r[istop] / initial_cond["Ms"]
    if Mcrat < -9.999:
        Mcrat = -9.999
    Lcrat = L_r[istop] / initial_cond["Ls"]
    if Lcrat < -9.999:
        Lcrat = -9.999

    solar_crat_table = [["","Solar", "Crat"],["Mass", initial_cond["Msolar"], Mcrat], ["Radius", initial_cond["Rsolar"], Rcrat], ["Luminosity", initial_cond["Lsolar"], Lcrat]]
    print(tabulate(solar_crat_table), headers="firstrow")
    mass_table = [["Hydrogen", "Helium", "Metals"], ["Mass Fractions", mass_fractions["X"], mass_fractions["Y"], mass_fractions["Z"]]]
    print(tabulate(mass_table), headers="firstrow")
    core_table = [["Rho", "Temperature", "Pressure", "Epsilon"], [rhocor, Tcore, Pcore, epscor]]
    print(tabulate(core_table, headers="firstrow"))

    # Print data from the cetner of the start outward, labeling convective or radiative zones zones.
    for ic in range(0 , istop):
        i = istop - ic + 1
        Qm = 1.0 - M_r[i]/initial_cond["Ms"]
        if dlPdlT[i] < constants["gamrat"]:
            rcf = 'c'
        else:
            rcf = 'r'
        if abs(dlPdlT[i]) > dlPlim:
            if dlPdlT[i] >= 0:
                dlPdlT = dlPlim
            else:
                dlPdlT = -dlPlim
            clim = '*'
        else:
            clim = ' '
        results = [["R [" + i + "]", r[i]], ["Qm", Qm], ["L_r[" + i + "]", L_r[i]], ["T[" + i + "", T[i]],
                   ["P["+ i + "]", P[i]], ["Rho[" + i + "]", rho[i]], ["Kappa[" + i + "]", kappa[i]],
                   ["Epsilon[" + i + "]", epslon[i]], ["dlPdlT[" + i + "]", dlPdlT[i]]]
        print(tabulate(results))
    print("**** The integration has been completed ****")

# Returns r, M_rip1, L_rip1, T_ip1, P_ip1
def STARTMDL(deltar, X, Z, mu, Rs, r_i, M_ri, L_ri, irc, returnArray):
    # r
    returnArray[0] = r_i + deltar
    # M_rip1
    returnArray[1] = M_ri
    # L_rip1
    returnArray[2] = L_ri

    if irc == 0:
        # T_ip1
        returnArray[3] = constants["G"] * returnArray[1]*mu*constants["m_H"]/(4.25*constants["k_B"])*(1/r - 1/Rs)
        A_bf = 4.34e+25 * Z * (1.0 + X)/constants["tog_bf"]
        A_ff = 3.68e+22 * constants["g_ff"] * (1.0 - Z) * (1.0 + X)
        Afac = A_bf + A_ff
        # L_ip1
        returnArray[4] = sqrt((1.0/4.25)*(16.0/3.0*pi*constants["a"]*constants["c"])*(constants["G"] * returnArray[1]/returnArray[2]) *
                     (constants["k_B"]/(Afac * mu * constants["m_H"]))) * returnArray[3]**4.25

    # This is the convective approximation
    else:
        returnArray[3] = constants["G"] * returnArray[1] * mu * constants["m_H"] / constants["k_B"] * (1/returnArray[0] - 1/Rs)/constants["gamrat"]
        returnArray[4] = constants["kPad"]*returnArray[3]**constants["gamrat"]

#calculates the values of density,opacity, guillotine-to-gaunt ration, energy gen rate
# for a radius r

# instead of passing in ierr going to use the return as an error code. 0 for success.
# Returns [rho, kappa, epslon]
def EOS(X, Z, XCNO, mu, P, T):

    oneo3 = .333333333333
    twoo3 = .666666666667

    # solve for density from the ideal gas law
    if T < 0 or P < 0:
       print('oops soupy')
       return 1
    Prad = constants['a']*(T**4)/3
    Pgas = P - Prad
    rho = (mu * constants['m_H']/constants['k_B']) * (Pgas/T)
    if rho < 0:
        print('oops soupy')
        return 1
    # Calc opacity, including guillatine-to-gaunt factor ratio
    tog_bf = 2.82 * (rho*(1+X))**.2
    k_bf = 4.34E25/tog_bf*Z(1 + X)* rho/(T**3.5)
    k_ff = 3.68E22*constants["g_ff"]*(1-Z)*(1+X)*rho/(T**3.5)
    k_e = .2*(1+X)
    kappa = k_bf + k_ff + k_e

    #calc eneregy generation by pp chain and CNO cycle
    #The screening factor for the pp chain is calculated as fpp

    T6 = T*1E-6
    fx = .133 * X * ((3+X)*rho)**.5 / (T6**1.5 + comp_eps)
    fpp = 1 + fx*X
    psipp = 1 + 1.412E8 * (1/x-1)*exp(-49.98*T6**-oneo3)
    Cpp = 1 + .0123*T6*oneo3 + .0109*T6**(-twoo3) - .000149*T6
    epspp = 2.38E6 * rho * x * X * fpp * psipp * Cpp * T6**(-twoo3) * exp(-33.8*T6**(-oneo3))
    CCNO = 1 + .0027*T6**oneo3 - .0078 * T6**twoo3 - .000149*T6
    epsCNO = 8.67E27 * rho * X * XCNO * CCNO * T6**(-twoo3) * exp(-152.28 * T6**(-oneo3))
    epslon = epspp + epsCNO

    return [rho, kappa, epslon]

# 'Hydrostatic equilibrium Pressure Gradient
def dPdr(r,M_r,rho):
    #dPdr
    return -constants["G"] * rho * M_r * r**-2

# 'Conservation of Mass
def dMdr(r,rho):
    #dMdr
    return 4 * pi * rho * r**2


# 'luminosity thingy
def dLdr(r, rho, epslon):
    #dLdr
    return 4 * pi * rho * epslon * r**2


# 'The temp one
def dTdr(r, M_r, L_r, T, rho, kappa, mu, irc):
    if irc == 0:
        #dTdr
        return - (3/(16 * pi * constants["a"] * constants["c"])) * kappa * rho * T**-3 * L_r * r**-2
    else:
        #dTdr
        return -1 / constants["gamrat"] * constants["G"] * M_r * r**-2 * mu * constants["m_H"] / constants["k_B"]

def RUNGE(f_im1, dfdr, f_i, r_im1, deltar, irc, X, Z, XCNO, mu, izone, ierr):
    f_temp = []
    df1 = []
    df2 = []
    df3 = []
    dr12 = deltar/2.0
    dr16 = deltar/6.0
    r12 = r_im1 + dr12
    r_i = r_im1 + deltar

    # Calculate intermediate derviatives from the fundamental stellar structure equations found in subroutine fundeq.
    for i in range(0,3):
        f_temp[i] = f_im1[i] + dr12*dfdr[i]
    FUNDEQ(r12, f_temp, df1, irc, X, Z, XCNO, mu, izone, ierr)
    for i in range(0, 3):
        f_temp[i] = f_im1[i] + dr12*df1[i]
    FUNDEQ(r12, f_temp, df2, irc, X, Z, XCNO, mu, izone, ierr)
    for i in range(0, 3):
        f_temp[i] = f_im1[i] + deltar*df2[i]
    FUNDEQ(r_i, f_temp, df3, irc, X, Z, XCNO, mu, izone, ierr)
    for i in range(0, 3):
        f_i[i] = f_im1[i] + dr16*(dfdr[i] + 2*df1[i] + 2.0 * df2[i] + df3[i])

def FUNDEQ(r, f, dfdr, irc, X, Z, XCNO, mu, izone, ierr):
    P = f[0]
    M_r = f[1]
    L_r = f[2]
    T = f[3]
    [rho, kappa, epslon] = EOS(X, Z, XCNO, mu, P, T, constants["tog_bf"], izone, ierr)
    dfdr[0] = dPdr(r, M_r, rho)
    dfdr[1] = dMdr(r, rho)
    dfdr[2] = dLdr(r, rho, epslon)
    dfdr[3] = dTdr(r, M_r, L_r, T, rho, kappa, mu, irc)
    
STATSTAR()