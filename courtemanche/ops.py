"""
ops.py — mathematical core of the model.

This module provides functions to compute the model equations,
as well as functions to retrieve default parameters and initial
values for the state variables.

This model describes the ionic currents and action potential dynamics of human atrial myocytes. 
It includes detailed formulations for major ionic currents (fast sodium current, L-type calcium current, inward rectifier potassium current, transient outward potassium current, 
rapid and slow delayed rectifier potassium currents, and Na⁺/Ca²⁺ exchanger), as well as calcium handling mechanisms.

The Courtemanche model is widely used as a reference atrial electrophysiology model. It has served as the basis for many subsequent atrial modeling studies, 
including investigations of atrial fibrillation and drug effects.

References:
Courtemanche M, Ramirez RJ, Nattel S. Ionic mechanisms underlying human atrial action potential properties: insights from a mathematical model. Am J Physiol. 1998 Jul;275(1):H301-21. 

DOI: 10.1152/ajpheart.1998.275.1.H301
"""

__all__ = (
    "get_variables",
    "get_parameters",
    "calc_rhs",  
    "calc_gating_variable",
    "calc_gating_variable_rush_larsen",
    "calc_cmdn",
    "calc_trpn",
    "calc_csqn",
    "calc_dnai",
    "calc_dki",
    "calc_dcai",
    "calc_dcaup",
    "calc_dcarel",
    "calc_equilibrum_potentials",
    "calc_ina",
    "calc_gating_m",
    "calc_gating_h",
    "calc_gating_j",
    "calc_ik1",
    "calc_ito",
    "calc_ikur",
    "calc_ikr",
    "calc_iks",
    "calc_ical",
    "calc_inak",
    "calc_inaca",
    "calc_ibca",
    "calc_ibna",
    "calc_ipca",
    "calc_irel",
    "calc_itr",
    "calc_iup",
    "calc_iupleak",
)

import math


def get_variables() -> dict[str, float]:
    """
    Returns default initial values for state variables.
    """
    return {
        "u": -81.2,
        "h": 0.965,
        "d": 0.000137,
        "xr": 0.0000329,
        "nai": 11.2,
        "ki": 139.0,
        "carel": 1.49,
        "oi": 0.999,
        "ui": 0.999,
        "vrel": 1.0,
        "m": 0.00291,
        "j": 0.978,
        "f": 0.999,
        "xs": 0.0187,
        "cai": 0.000102,
        "caup": 1.49,
        "oa": 0.0304,
        "ua": 0.00496,
        "fca": 0.775,
        "irel": 0.0,
        "urel": 0.0,
        "wrel": 0.999,
    }


def get_parameters() -> dict[str, float]:
    """
    Returns default parameter values for the model.
    """
    return {
        "R": 8.3143,
        "T": 310.0,
        "F": 96.4867,
        "Cm": 100.0,
        "Vc": 20100.0,
        "Vj": 13668.0,
        "Vup": 1109.52,
        "Vrel": 96.48,
        "ko": 5.4,
        "nao": 140.0,
        "cao": 1.8,
        "gna": 7.8,
        "gk1": 0.09,
        "gto": 0.1652,
        "gkr": 0.0294,
        "gks": 0.129,
        "gcal": 0.1238,
        "gcab": 0.00113,
        "gnab": 0.000674,
        "inakmax": 0.6,
        "inacamax": 1600.0,
        "ipcamax": 0.275,
        "iupmax": 0.005,
        "kq10": 3.0,
        "gamma": 0.35, # for I_NaCa
        "kmnai": 10.0,
        "kmko": 1.5,
        "kmnancx": 87.5,
        "kmcancx": 1.38,
        "ksatncx": 0.1,
        "krel": 30.0,
        "kup": 0.00092,
        "caupmax": 15.0,
        "cmdnmax": 0.05,
        "trpnmax": 0.07,
        "csqnmax": 10.0,
        "kmcmdn": 0.00238,
        "kmtrpn": 0.0005,
        "kmcsqn": 0.8,        
        "ibk": 0.0,
    }



def calc_rhs(ina, ik1, ito, ikur, ikr, iks, ical, ipca, inak, inaca, ibna, ibca, Cm):
    """
    Computes the ionic currents density (in pA/pF) for the model.

    Parameters
    ----------
    ina : float
        Fast sodium current.
    ik1 : float
        Time-independent potassium current.
    ito : float
        Transient outward potassium current
    ikur : float
        Ultra-rapid delayed rectifier potassium current.
    ikr : float
        Rapid delayed rectifier potassium current.
    iks : float
        Slow delayed rectifier potassium current.
    ical : float
        L-type calcium current.
    ipca : float
        Sarcolemmal calcium pump current.
    inak : float
        Sodium-potassium pump current.
    inaca : float
        Sodium-calcium exchanger current.
    ibna : float
        Background sodium current.
    ibca : float
        Background calcium current.
    Cm : float
        Cell membrane capacitance.
    """
    rhs =  ina + ik1 + ito + ikur + ikr + iks + ical + ipca + inak + inaca + ibna + ibca
    return rhs / Cm

def calc_gating_variable(x, x_inf, tau_x,):
    """
    Calculates the gating variable using the steady-state value and time constant.

    Parameters
    ----------
    x : float
        Current value of the gating variable.
    x_inf : float
        Steady-state value of the gating variable.
    tau_x : float
        Time constant for the gating variable (ms).
    """
    return (x_inf - x) / tau_x

def calc_gating_variable_rush_larsen(x, x_inf, tau_x, dt, exp=math.exp):
    """
    Calculates the gating variable using the Rush-Larsen method.

    Parameters
    ----------
    x : float
        Current value of the gating variable.
    x_inf : float
        Steady-state value of the gating variable.
    tau_x : float
        Time constant for the gating variable (ms).
    exp : callable
        Exponential function to use (default: math.exp).
    """
    return x_inf - (x_inf - x)*exp(-dt/tau_x)


def calc_cmdn(cmdnmax, kmcmdn, cai):
    """
    Calculates the concentration of calmodulin.

    Parameters
    ----------
    cmdnmax : float
        Maximum concentration of calmodulin.
    kmcmdn : float
        Dissociation constant for calmodulin.
    cai : float
        Intracellular calcium concentration.
    """
    cmdn = cmdnmax*cai/(cai + kmcmdn)
    return cmdn


def calc_trpn(trpnmax, kmtrpn, cai):
    """
    Calculates the concentration of troponin.

    Parameters
    ----------
    trpnmax : float
        Maximum concentration of troponin.
    kmtrpn : float
        Dissociation constant for troponin.
    cai : float
        Intracellular calcium concentration.
    """
    trpn = trpnmax*cai/(cai + kmtrpn)
    return trpn


def calc_csqn(csqnmax, kmcsqn, carel):
    """
    Calculates the concentration of calsequestrin.

    Parameters
    ----------
    csqnmax : float
        Maximum concentration of calsequestrin.
    kmcsqn : float
        Dissociation constant for calsequestrin.
    carel : float
        Calcium concentration in the release compartment.
    """
    csqn = csqnmax*carel/(carel + kmcsqn)
    return csqn


def calc_dnai(inak, inaca, ibna, ina, F, Vj):
    """
    Calculates the intracellular sodium concentration.

    Parameters
    ----------
    inak : float
        Sodium-potassium pump current.
    inaca : float
        Sodium-calcium exchanger current.
    ibna : float
        Background sodium current.
    ina : float
        Fast sodium current.
    F : float
        Faraday's constant.
    Vj : float
        Cell volume.
    """
    dnai = (-3*inak-3*inaca - ibna - ina)/(F*Vj)
    return dnai


def calc_dki(inak, ik1, ito, ikur, ikr, iks, ibk, F, Vj):

    """
    Calculates the intracellular potassium concentration.

    Parameters
    ----------
    inak : float
        Sodium-potassium pump current.
    ik1 : float
        Time-independent potassium current.
    ito : float
        Transient outward potassium current.
    ikur : float
        Ultra-rapid delayed rectifier potassium current.
    ikr : float
        Rapid delayed rectifier potassium current.
    iks : float
        Slow delayed rectifier potassium current.
    ibk : float
        Background potassium current.
    F : float
        Faraday's constant.
    Vj : float
        Cell volume.
    """
    dki = (2*inak - ik1 - ito - ikur - ikr - iks - ibk)/(F*Vj)
    return dki


def calc_dcai(cai, inaca, ipca, ical, ibca, iup, iupleak, irel, Vrel, Vup, trpnmax, kmtrpn, cmdnmax, kmcmdn, F, Vj): 
    """
    Calculates the intracellular calcium concentration.

    Parameters
    ----------
    cai : float
        Intracellular calcium concentration.
    inaca : float
        Sodium-calcium exchanger current.
    ipca : float
        Sarcolemmal calcium pump current.
    ical : float
        L-type calcium current.
    ibca : float
        Background calcium current.
    iup : float
        Uptake of calcium into the NSR.
    iupleak : float
        Leak of calcium from the NSR.
    irel : float
        Calcium release from the JSR.
    Vrel : float
        Volume of the release compartment.
    Vup : float
        Volume of the uptake compartment.
    trpnmax : float
        Maximum concentration of troponin.
    kmtrpn : float
        Dissociation constant for troponin.
    cmdnmax : float
        Maximum concentration of calmodulin.
    kmcmdn : float
        Dissociation constant for calmodulin.
    F : float
        Faraday's constant.
    Vj : float
        Cell volume.
    """
    B1 = (2*inaca - ipca - ical - ibca)/(2*F*Vj) + (Vup*(iupleak - iup) + irel*Vrel)/Vj
    B2 = 1 + (trpnmax*kmtrpn)/((cai + kmtrpn)**2) + (cmdnmax*kmcmdn)/((cai + kmcmdn)**2)
    dcai = B1/B2
    return dcai


def calc_dcaup(iup, iupleak, itr, Vrel, Vup):
    """
    Calculates the calcium concentration in the up compartment.

    Parameters
    ----------
    iup : float
        Uptake of calcium into the NSR.
    iupleak : float
        Leak of calcium from the NSR.
    itr : float
        Transfer of calcium from the NSR to the JSR.
    Vrel : float
        Volume of the release compartment.
    Vup : float
        Volume of the uptake compartment.
    """
    dcaup = iup - iupleak - itr*(Vrel/Vup)
    return dcaup


def calc_dcarel(carel, itr, irel, csqnmax, kmcsqn):
    """
    Calculates the calcium concentration in the release compartment.

    Parameters
    ----------
    carel : float
        Calcium concentration in the release compartment.
    itr : float
        Transfer of calcium from the NSR to the JSR.
    irel : float
        Calcium release from the JSR.
    csqnmax : float
        Maximum concentration of calsequestrin.
    kmcsqn : float
        Dissociation constant for calsequestrin.
    """
    dcarel = (itr - irel)/(1 + (csqnmax*kmcsqn)/((carel + kmcsqn)**2))
    return dcarel


def calc_equilibrum_potentials(nai, nao, ki, ko, cai, cao, R, T, F, log=math.log):
    """
    Calculates the equilibrum potentials for the cell.

    Parameters
    ----------
    nai : float
        Intracellular sodium concentration.
    nao : float
        Extracellular sodium concentration.
    ki : float
        Intracellular potassium concentration.
    ko : float
        Extracellular potassium concentration.
    cai : float
        Intracellular calcium concentration.
    cao : float
        Extracellular calcium concentration.
    R : float
        Universal gas constant.
    T : float
        Absolute temperature.
    F : float
        Faraday's constant.
    log : callable
        Logarithm function to use (default: math.log).
    """
    ena = (R*T/F)*log(nao/nai)

    ek = (R*T/F)*log(ko/ki)

    safe_cai = max(cai, 1e-7)
    
    eca = (R*T/(2*F))*log(cao/safe_cai)
    return ena, ek, eca


def calc_ina(u, m, h, j, gna, ena, Cm):
    """
    Calculates the fast sodium current.

    Parameters
    ----------
    u : float
        Membrane potential.
    m : float
        Gating variable m.
    h : float
        Gating variable h.
    j : float
        Gating variable j.
    gna : float
        Maximum conductance for the fast sodium current.
    ena : float
        Equilibrium potential for sodium.
    Cm : float
        Cell membrane capacitance.
    """
    ina = Cm *  gna*(m**3)*h*j*(u - ena)
    return ina


def calc_gating_m(m, u, dt, exp=math.exp):
    """
    Calculates the gating variable m for the fast sodium current.

    Parameters
    ----------
    m : float
        Current value of the gating variable m.
    u : float
        Membrane potential.
    dt : float
        Time step (ms).
    exp : callable
        Exponential function to use (default: math.exp).
    """
    am = 0
    if u == -47.13:
        am = 3.2
    else:
        am = 0.32 * (u + 47.13) / (1 - exp(-0.1 * (u + 47.13)))

    bm = 0.08*exp(-u/11)
    m_inf = am/(am + bm)
    tau_m = 1/(am + bm)
    m = calc_gating_variable_rush_larsen(m, m_inf, tau_m, dt)

    return m


def calc_gating_h(h, u, dt, exp=math.exp):
    """
    Calculates the gating variable h for the fast sodium current.

    Parameters
    ----------
    h : float
        Current value of the gating variable h.
    u : float
        Membrane potential.
    dt : float
        Time step (ms).
    exp : callable
        Exponential function to use (default: math.exp).
    """
    ah = 0.135*exp(-(80 + u)/6.8)
    bh = 3.56*exp(0.079*u) + 310000*exp(0.35*u)
    if u >= -40:
        ah = 0    
        bh = 1/(0.13*(1 + exp(-(u + 10.66)/11.1)))
    
    h_inf = ah/(ah + bh)
    tau_h = 1/(ah + bh)
    h = calc_gating_variable_rush_larsen(h, h_inf, tau_h, dt)
    return h


def calc_gating_j(j, u, dt, exp=math.exp):
    """
    Calculates the gating variable j for the fast sodium current.

    Parameters
    ----------
    j : float
        Current value of the gating variable j.
    u : float
        Membrane potential.
    dt : float
        Time step (ms).
    exp : callable
        Exponential function to use (default: math.exp).
    """
    
    aj = (-127140*exp(0.2444*u) - 0.00003474*exp(-0.04391*u))*(u + 37.78)/(1 + exp(0.311*(u + 79.23)))
    bj = 0.1212*exp(-0.01052*u)/(1 + exp(-0.1378*(u + 40.14)))
    if u >= -40:
        aj = 0
        bj = 0.3*exp(-0.0000002535*u)/(1 + exp(-0.1*(u + 32)))
    j_inf = aj/(aj + bj)
    tau_j = 1/(aj + bj)
    j = calc_gating_variable_rush_larsen(j, j_inf, tau_j, dt)
    return j


def calc_ik1(u, gk1, ek, Cm, exp=math.exp):
    """
    Calculates the time-independent potassium current.

    Parameters
    ----------
    u : float
        Membrane potential.
    gk1 : float
        Maximum conductance for the time-independent potassium current.
    ek : float
        Equilibrium potential for potassium.
    Cm : float
        Cell membrane capacitance.
    exp : callable
        Exponential function to use (default: math.exp).
    """
    ik1 = Cm *  gk1*(u - ek)/(1 + exp(0.07*(u + 80)))
    return ik1


def calc_ito(u, dt, kq10, oa, oi, gto, ek, Cm, exp=math.exp):
    """
    Calculates the transient outward potassium current.

    Parameters
    ----------
    u : float
        Membrane potential.
    dt : float
        Time step (ms).
    kq10 : float
        Temperature adjustment factor.
    oa : float
        Gating variable oa.
    oi : float
        Gating variable oi.
    gto : float
        Maximum conductance for the transient outward potassium current.
    ek : float
        Equilibrium potential for potassium.
    Cm : float
        Cell membrane capacitance.
    exp : callable
        Exponential function to use (default: math.exp).
    """
    ao = 0.65/(exp(-(u + 10)/8.5) + exp(-(u - 30)/59.0))
    bo = 0.65/(2.5 + exp((u + 82)/17.0))

    tau_o = 1/(kq10*(ao + bo))
    o_inf = 1/(1 + exp(-(u + 20.47)/17.54))

    aoi = 1/(18.53 + exp((u + 113.7)/10.95))
    boi = 1/(35.56 + exp(-(u + 1.26)/7.44))

    tau_oi = 1/(kq10*(aoi + boi))
    oi_inf = 1/(1 + exp((u + 43.1)/5.3))

    oa = calc_gating_variable_rush_larsen(oa, o_inf, tau_o, dt)
    oi = calc_gating_variable_rush_larsen(oi, oi_inf, tau_oi, dt)

    ito = Cm *  gto*(oa**3)*oi*(u - ek)  

    return ito, oa, oi


def calc_ikur(u, dt, kq10, ua, ui, ek, Cm, exp=math.exp):
    """
    Calculates the ultra-rapid delayed rectifier potassium current.

    Parameters
    ----------
    u : float
        Membrane potential.
    dt : float
        Time step (ms).
    kq10 : float
        Temperature adjustment factor.
    ua : float
        Gating variable ua.
    ui : float
        Gating variable ui.
    ek : float
        Equilibrium potential for potassium.
    Cm : float
        Cell membrane capacitance.
    exp : callable
        Exponential function to use (default: math.exp).
    """
    gkur = 0.005 + 0.05/(1 + exp(-(u - 15)/13.0))

    aua = 0.65/(exp(-(u + 10)/8.5) + exp(-(u - 30)/59.0))
    bua = 0.65/(2.5 + exp((u + 82)/17.0))
    tau_ua = 1/(kq10*(aua + bua))
    ua_inf = 1/(1 + exp(-(u + 30.3)/9.6))
    aui = 1/(21 + exp(-(u - 185)/28.0))
    bui = exp((u - 158)/16.0)

    tau_ui = 1/(kq10*(aui + bui))
    ui_inf = 1/(1 + exp((u - 99.45)/27.48))

    ua = calc_gating_variable_rush_larsen(ua, ua_inf, tau_ua, dt)
    ui = calc_gating_variable_rush_larsen(ui, ui_inf, tau_ui, dt)

    ikur = Cm *  gkur*(ua**3)*ui*(u - ek)

    return ikur, ua, ui


def calc_ikr(u, dt, xr, gkr, ek, Cm, exp=math.exp):
    """
    Calculates the rapid delayed rectifier potassium current.

    Parameters
    ----------
    u : float
        Membrane potential.
    dt : float
    xr : float
        Gating variable xr.
    gkr : float
        Maximum conductance for the rapid delayed rectifier potassium current.
    ek : float
        Equilibrium potential for potassium.
    Cm : float
        Cell membrane capacitance.
    exp : callable
        Exponential function to use (default: math.exp).
    """
    gkr = 0.0294 # * np.sqrt(ko / 5.4)
    axr = 0.0003*(u + 14.1)/(1 - exp(-(u + 14.1)/5))
    bxr = 0.000073898*(u - 3.3328)/(exp((u - 3.3328)/5.1237) - 1)

    tau_xr = 1/(axr + bxr)
    xr_inf = 1/(1 + exp(-(u + 14.1)/6.5))

    xr = calc_gating_variable_rush_larsen(xr, xr_inf, tau_xr, dt)

    ikr = Cm * (gkr*xr*(u - ek))/(1 + exp((u + 15)/22.4))

    return ikr, xr


def calc_iks(u, dt, xs, gks, ek, Cm, exp=math.exp, sqrt=math.sqrt):
    """
    Calculates the slow delayed rectifier potassium current.
    """
    axs = 0.00004*(u - 19.9)/(1 - exp(-(u - 19.9)/17))
    bxs = 0.000035*(u - 19.9)/(exp((u - 19.9)/9) - 1)

    tau_xs = 1/(2*(axs + bxs))
    xs_inf = 1/sqrt(1 + exp(-(u - 19.9)/12.7))

    xs = calc_gating_variable_rush_larsen(xs, xs_inf, tau_xs, dt)

    iks = Cm *  gks*(xs**2)*(u - ek)

    return iks, xs


def calc_ical(u, dt, d, f, cai, gcal, fca, Cm, exp=math.exp):
    """
    Calculates the L-type calcium current.

    Parameters
    ----------
    u : float
        Membrane potential.
    dt : float
    d : float
        Gating variable d.
    f : float
        Gating variable f.
    cai : float
        Intracellular calcium concentration.
    gcal : float
        Maximum conductance for the L-type calcium current.
    fca : float
        Gating variable fca.
    Cm : float
        Cell membrane capacitance.      
    exp : callable
        Exponential function to use (default: math.exp).
    """
    tau_d = (1 - exp(-(u + 10)/6.24))/(0.035*(u + 10)*(1 + exp(-(u + 10)/6.24)))
    d_inf = 1/(1 + exp(-(u + 10)/8.0))

    tau_f = 9/(0.0197*exp(-(0.0337**2)*((u + 10)**2)) + 0.02)
    f_inf = 1/(1 + exp((u + 28)/6.9))

    tau_fca = 2
    fca_inf = 1/(1 + cai/0.00035)

    d   = calc_gating_variable_rush_larsen(d, d_inf, tau_d, dt)
    f   = calc_gating_variable_rush_larsen(f, f_inf, tau_f, dt)
    fca = calc_gating_variable_rush_larsen(fca, fca_inf, tau_fca, dt)

    ical = Cm *  gcal*d*f*fca*(u - 65) 

    return ical, d, f, fca


def calc_inak(inakmax, nai, nao, ko, kmnai, kmko, F, u, R, T, Cm, exp=math.exp):
    """
    Calculates the sodium-potassium pump current.

    Parameters
    ----------
    inakmax : float
        Maximum current for the sodium-potassium pump.
    nai : float
        Intracellular sodium concentration.
    nao : float
        Extracellular sodium concentration.
    ko : float
        Extracellular potassium concentration.
    kmnai : float
        Dissociation constant for intracellular sodium.
    kmko : float
        Dissociation constant for extracellular potassium.
    F : float
        Faraday's constant.
    u : float
        Membrane potential.
    R : float
        Universal gas constant.
    T : float
        Absolute temperature.
    Cm : float
        Cell membrane capacitance.
    exp : callable
        Exponential function to use (default: math.exp).
    """
    s = (1/7.0)*(exp(nao/67.3) - 1)
    fnak = 1/(1 + 0.1245*exp(-0.1*(F*u)/(R*T)) + 0.0365*s*exp(-(F*u)/(R*T)))

    inak = Cm *  inakmax*fnak*(1/(1 + (kmnai/nai)**1.5))*(ko/(ko + kmko))

    return inak


def calc_inaca(inacamax, nai, nao, cai, cao, kmnancx, kmcancx, ksatncx, F, u, R, T, Cm, exp=math.exp):
    """
    Calculates the sodium-calcium exchanger current.

    Parameters
    ----------
    inacamax : float
        Maximum current for the sodium-calcium exchanger.
    nai : float
        Intracellular sodium concentration.
    nao : float
        Extracellular sodium concentration.
    cai : float
        Intracellular calcium concentration.
    cao : float
        Extracellular calcium concentration.
    kmnancx : float
        Dissociation constant for intracellular sodium.
    kmcancx : float
        Dissociation constant for extracellular calcium.
    ksatncx : float
        Saturation factor for the sodium-calcium exchanger.
    F : float
        Faraday's constant.
    u : float
        Membrane potential.
    R : float
        Universal gas constant.
    T : float
        Absolute temperature.
    Cm : float
        Cell membrane capacitance.
    exp : callable
        Exponential function to use (default: math.exp).
    """
    gamma = 0.35

    # Exponential terms with clamping
    exp_term = exp(gamma * (F * u) / (R * T))
    exp_rev_term = exp((gamma - 1) * (F * u) / (R * T))

    # Numerator
    numerator = inacamax * (exp_term * nai**3 * cao - exp_rev_term * nao**3 * cai)

    # Denominator
    term1 = (kmnancx**3 + nao**3)  # (K_m,Na^3 + [Na+]_i^3)
    term2 = (kmcancx + cao)        # (K_m,Ca + [Ca2+]_o)
    term3 = (1 + ksatncx * exp_rev_term)  # (1 + k_sat * exp(...))

    denominator = term1 * term2 * term3

    # Calculate INaCa
    inaca = Cm *  numerator / denominator

    return inaca


def calc_ibca(gcab, eca, u, Cm):
    """
    Calculates the background calcium current.

    Parameters
    ----------
    gcab : float
        Maximum conductance for the background calcium current.
    eca : float
        Equilibrium potential for calcium.
    u : float
        Membrane potential.
    Cm : float
        Cell membrane capacitance.
    """
    ibca = Cm *  gcab*(u - eca)
    return ibca


def calc_ibna(gnab, ena, u, Cm):
    """
    Calculates the background sodium current.

    Parameters
    ----------
    gnab : float
        Maximum conductance for the background sodium current.
    ena : float
        Equilibrium potential for sodium.
    u : float
        Membrane potential.
    Cm : float
        Cell membrane capacitance.
    """
    ibna = Cm *  gnab*(u - ena)
    return ibna

def calc_ipca(ipcamax, cai, Cm):
    """
    Calculates the sarcolemmal calcium pump current.

    Parameters
    ----------
    ipcamax : float
        Maximum current for the sarcolemmal calcium pump.
    cai : float
        Intracellular calcium concentration.
    Cm : float
        Cell membrane capacitance.
    """
    ipca = Cm *  ipcamax*cai/(cai + 0.0005)
    return ipca

def calc_irel(dt, urel, vrel, irel, wrel, ical, inaca, krel, carel, cai, u, F, Vrel, exp=math.exp): 
    """
    Calculates the calcium release from the JSR.

    Parameters
    ----------
    dt : float
        Time step (ms).
    urel : float    
        Gating variable urel.
    vrel : float
        Gating variable vrel.
    irel : float
        Calcium release from the JSR.
    wrel : float
        Gating variable wrel.
    ical : float
        L-type calcium current.
    inaca : float
        Sodium-calcium exchanger current.
    krel : float
        Scaling factor for the calcium release.
    carel : float
        Calcium concentration in the release compartment.
    cai : float
        Intracellular calcium concentration.
    u : float
        Membrane potential.
    F : float
        Faraday's constant.
    Vrel : float
        Volume of the release compartment.
    exp : callable
        Exponential function to use (default: math.exp).
    """
    tau_u = 8

    Fn = 1e-12*Vrel*irel - ((5*1e-13)/F)*(0.5*ical - 0.2*inaca) 

    u_inf = 1/(1 + exp(-(Fn - 3.4175e-13)/13.67e-16))

    tau_v = 1.91 + 2.09/(1 + exp(-(Fn - 3.4175e-13)/13.67e-16))
    v_inf = 1 - 1/(1 + exp(-(Fn - 6.835e-14)/13.67e-16))

    tau_w = 6 * (1 - exp(-(u - 7.9) / 5.0)) / ((1 + 0.3 * exp(-(u - 7.9) / 5.0)) * (u - 7.9))
    w_inf = 1 - 1/(1 + exp(-(u - 40)/17.0))

    urel = calc_gating_variable_rush_larsen(urel, u_inf, tau_u, dt)
    vrel = calc_gating_variable_rush_larsen(vrel, v_inf, tau_v, dt)
    wrel = calc_gating_variable_rush_larsen(wrel, w_inf, tau_w, dt)

    irel = krel*(urel**2)*vrel*wrel*(carel - cai)

    return irel, urel, vrel, wrel

def calc_itr(caup, carel):
    """
    Calculates the transfer of calcium from the NSR to the JSR.

    Parameters
    ----------
    caup : float
        Calcium concentration in the uptake compartment.
    carel : float
        Calcium concentration in the release compartment.
    """
    tautr = 180
    itr = (caup - carel)/tautr
    return itr

def calc_iup(iupmax, cai, kup):
    """
    Calculates the uptake of calcium into the NSR.

    Parameters
    ----------
    iupmax : float
        Maximum uptake of calcium into the NSR.
    cai : float
        Intracellular calcium concentration.
    kup : float
        Dissociation constant for the uptake.
    """
    iup = iupmax/(1 + (kup/cai))
    return iup

def calc_iupleak(caup, caupmax, iupmax):
    """
    Calculates the leak of calcium from the NSR.

    Parameters
    ----------
    caup : float
        Calcium concentration in the uptake compartment.
    caupmax : float
        Maximum calcium concentration in the uptake compartment.
    iupmax : float
        Maximum uptake of calcium into the NSR.
    """
    iupleak = (caup/caupmax)*iupmax
    return iupleak