"""
py — mathematical core of the model.

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
    "calc_where",
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
    "calc_ena",
    "calc_ek",
    "calc_eca",
    "calc_ina",
    "calc_am",
    "calc_bm",
    "calc_tau",
    "calc_inf",
    "calc_ah",
    "calc_bh",
    "calc_aj",
    "calc_bj",
    "calc_ik1",
    "calc_tau_oa",
    "calc_oa_inf",
    "calc_tau_oi",
    "calc_oi_inf",
    "calc_ito",
    "calc_tau_ua",
    "calc_ua_inf",
    "calc_tau_ui",
    "calc_ui_inf",
    "calc_ikur",
    "calc_tau_xr",
    "calc_xr_inf",
    "calc_ikr",
    "calc_tau_xs",
    "calc_xs_inf",
    "calc_iks",
    "calc_tau_d",
    "calc_d_inf",
    "calc_tau_f",
    "calc_f_inf",
    "calc_tau_fca",
    "calc_fca_inf",
    "calc_ical",
    "calc_inak",
    "calc_inaca",
    "calc_ibca",
    "calc_ibna",
    "calc_ipca",
    "calc_Fn",
    "calc_tau_urel",
    "calc_urel_inf",
    "calc_tau_vrel",
    "calc_vrel_inf",
    "calc_tau_wrel",
    "calc_wrel_inf",
    "calc_irel",
    "calc_itr",
    "calc_iup",
    "calc_iupleak",
    
)

from math import exp, log, sqrt 


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
        # "irel": 0.0,
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


def ionic_step(dt, u, m, h, j, oa, oi, ua, ui, xr, xs, d, f, fca, urel, vrel,
               wrel, caup, carel, cai, nai, ki,
               R, T, F, Cm, Vc, Vj, Vup, Vrel, ko, nao, cao, gna, gk1, gto,
               gkr, gks, gcal, gcab, gnab, inakmax, inacamax, ipcamax, iupmax,
               kq10, gamma, kmnai, kmko, kmnancx, kmcancx, ksatncx, krel, kup,
               caupmax, cmdnmax, trpnmax, csqnmax, kmcmdn, kmtrpn, kmcsqn, ibk,):
    
    """Performs a single time step update for the model state variables.
    
    Parameters
    ----------
    dt : float
        Time step (ms).
    u : float
        Membrane potential (mV).
    m, h, j : float
        Gating variables for the fast sodium current.
    oa, oi : float
        Gating variables for the transient outward potassium current.
    ua, ui : float
        Gating variables for the ultra-rapid delayed rectifier potassium current.
    xr : float
        Gating variable for the rapid delayed rectifier potassium current.
    xs : float
        Gating variable for the slow delayed rectifier potassium current.
    d, f : float
        Gating variable for the L-type calcium current.
    fca : float
        Gating variable for the calcium release current.
    urel, vrel, wrel : float
        Gating variables for the calcium release current.
    caup : float
        Calcium concentration in the sarcoplasmic reticulum.
    carel : float
        Calcium concentration in the cytosol.
    cai : float
        Intracellular calcium concentration.
    nai : float
        Intracellular sodium concentration.
    ki : float
        Intracellular potassium concentration.
    """

    ena = calc_ena(nai, nao, R, T, F)
    ek = calc_ek(ki, ko, R, T, F)
    eca = calc_eca(cai, cao, R, T, F)

    am = calc_am(u)
    bm = calc_bm(u)
    tau_m = calc_tau(am, bm)
    m_inf = calc_inf(am, bm)
    m_new = calc_gating_variable_rush_larsen(m, m_inf, tau_m, dt)

    ah = calc_ah(u)
    bh = calc_bh(u) 
    tau_h = calc_tau(ah, bh)
    h_inf = calc_inf(ah, bh)
    h_new = calc_gating_variable_rush_larsen(h, h_inf, tau_h, dt)

    aj = calc_aj(u)
    bj = calc_bj(u)
    tau_j = calc_tau(aj, bj)
    j_inf = calc_inf(aj, bj)
    j_new = calc_gating_variable_rush_larsen(j, j_inf, tau_j, dt)

    ina = calc_ina(u, m, h, j, gna, ena, Cm)

    ik1 = calc_ik1(u, gk1, ek, Cm)

    tau_oa = calc_tau_oa(u, kq10)
    oa_inf = calc_oa_inf(u)
    oa_new = calc_gating_variable_rush_larsen(oa, oa_inf, tau_oa, dt)

    tau_oi = calc_tau_oi(u, kq10)
    oi_inf = calc_oi_inf(u)
    oi_new = calc_gating_variable_rush_larsen(oi, oi_inf, tau_oi, dt)

    ito = calc_ito(u, oa, oi, gto, ek, Cm)

    tau_ua = calc_tau_ua(u, kq10)
    ua_inf = calc_ua_inf(u)
    ua_new = calc_gating_variable_rush_larsen(ua, ua_inf, tau_ua, dt)

    tau_ui = calc_tau_ui(u, kq10)
    ui_inf = calc_ui_inf(u)
    ui_new = calc_gating_variable_rush_larsen(ui, ui_inf, tau_ui, dt)

    ikur = calc_ikur(u, ua, ui, ek, Cm)

    tau_xr = calc_tau_xr(u)
    xr_inf = calc_xr_inf(u)
    xr_new = calc_gating_variable_rush_larsen(xr, xr_inf, tau_xr, dt)

    ikr = calc_ikr(u, xr, gkr, ek, Cm)

    tau_xs = calc_tau_xs(u)
    xs_inf = calc_xs_inf(u)
    xs_new = calc_gating_variable_rush_larsen(xs, xs_inf, tau_xs, dt)

    iks = calc_iks(u, xs, gks, ek, Cm)

    tau_d = calc_tau_d(u)
    d_inf = calc_d_inf(u)
    d_new = calc_gating_variable_rush_larsen(d, d_inf, tau_d, dt)

    tau_f = calc_tau_f(u)
    f_inf = calc_f_inf(u)
    f_new = calc_gating_variable_rush_larsen(f, f_inf, tau_f, dt)

    tau_fca = calc_tau_fca()
    fca_inf = calc_fca_inf(cai)
    fca_new = calc_gating_variable_rush_larsen(fca, fca_inf, tau_fca, dt)

    ical = calc_ical(u, d, f, gcal, fca, Cm)

    inak = calc_inak(inakmax, nai, nao, ko, kmnai, kmko, F, u, R, T, Cm)

    inaca = calc_inaca(inacamax, nai, nao, cai, cao, kmnancx, kmcancx, ksatncx, F, u, R, T, Cm)

    ibca = calc_ibca(gcab, eca, u, Cm)

    ibna = calc_ibna(gnab, ena, u, Cm)

    ipca = calc_ipca(ipcamax, cai, Cm)

    irel = calc_irel(urel, vrel, wrel, krel, carel, cai)

    Fn = calc_Fn(irel, ical, inaca, F, Vrel)

    tau_urel = calc_tau_urel()
    urel_inf = calc_urel_inf(Fn)
    urel_new = calc_gating_variable_rush_larsen(urel, urel_inf, tau_urel, dt)

    tau_vrel = calc_tau_vrel(Fn)
    vrel_inf = calc_vrel_inf(Fn)
    vrel_new = calc_gating_variable_rush_larsen(vrel, vrel_inf, tau_vrel, dt)

    tau_wrel = calc_tau_wrel(u)
    wrel_inf = calc_wrel_inf(u)
    wrel_new = calc_gating_variable_rush_larsen(wrel, wrel_inf, tau_wrel, dt)


    itr = calc_itr(caup, carel)
    iup = calc_iup(iupmax, cai, kup)
    iupleak = calc_iupleak(caup, caupmax, iupmax)

    dcaup = calc_dcaup(iup, iupleak, itr, Vrel, Vup)

    dnai = calc_dnai(inak, inaca, ibna, ina, F, Vj)

    dki = calc_dki(inak, ik1, ito, ikur, ikr, iks, ibk, F, Vj)

    dcai = calc_dcai(cai, inaca, ipca, ical, ibca, iup, iupleak, irel, Vrel, Vup,
                     trpnmax, kmtrpn, cmdnmax, kmcmdn, F, Vj)

    dcarel = calc_dcarel(carel, itr, irel, csqnmax, kmcsqn)

    rhs = - calc_rhs(ina, ik1, ito, ikur, ikr, iks, ical, ipca, inak, inaca, ibna, ibca, Cm)

    caup_new = caup + dt * dcaup
    nai_new = nai + dt * dnai
    ki_new = ki + dt * dki
    cai_new = cai + dt * dcai
    carel_new = carel + dt * dcarel

    return (rhs, m_new, h_new, j_new, oa_new, oi_new, ua_new, ui_new, xr_new, xs_new,
            d_new, f_new, fca_new, urel_new, vrel_new, wrel_new,
            caup_new, carel_new, cai_new, nai_new, ki_new)


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


def calc_where(cond, x, y):
    if cond:
        return x
    return y


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

def calc_gating_variable_rush_larsen(x, x_inf, tau_x, dt):
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


def calc_ena(nai, nao, R, T, F):
    """
    Calculates the equilibrium potential for sodium.

    Parameters
    ----------
    nai : float
        Intracellular sodium concentration.
    nao : float
        Extracellular sodium concentration.
    R : float
        Universal gas constant.
    T : float
        Absolute temperature.
    F : float
        Faraday's constant.
    """
    ena = (R*T/F)*log(nao/nai)
    return ena


def calc_ek(ki, ko, R, T, F):
    """
    Calculates the equilibrium potential for potassium.

    Parameters
    ----------
    ki : float
        Intracellular potassium concentration.
    ko : float
        Extracellular potassium concentration.
    R : float
        Universal gas constant.
    T : float
        Absolute temperature.
    F : float
        Faraday's constant.
    """
    ek = (R*T/F)*log(ko/ki)
    return ek


def calc_eca(cai, cao, R, T, F):
    """
    Calculates the equilibrium potential for calcium.

    Parameters
    ----------
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
    """
    # safe_cai = max(cai, 1e-7)
    eca = (R*T/(2*F))*log(cao/cai)
    return eca


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

    ina = Cm * gna*(m**3)*h*j*(u - ena)
    return ina


def calc_am(u):
    """
    Calculates the alpha rate for the gating variable m.

    Parameters
    ----------
    u : float
        Membrane potential.

    Note
    ----
        The singularity at u = -47.13 mV is handled using
        lim x->0 : x / (1 - exp(x)) = -1
    """

    am = calc_where(abs(u + 47.13) < 1e-5,
                    0.32 * 10,
                    0.32 * (u + 47.13) / (1 - exp(-0.1 * (u + 47.13))))
    return am


def calc_bm(u):
    """
    Calculates the beta rate for the gating variable m.

    Parameters
    ----------
    u : float
        Membrane potential.
    """

    bm = 0.08 * exp(-u / 11)
    return bm

def calc_tau(a, b):
    """
    Calculates the time constant for the gating variable.

    Parameters
    ----------
    a : float
        Alpha rate for the gating variable.
    b : float
        Beta rate for the gating variable.
    """

    tau = 1/(a + b)
    return tau


def calc_inf(a, b):
    """
    Calculates the steady-state value for the gating variable.

    Parameters
    ----------
    a : float
        Alpha rate for the gating variable.
    b : float
        Beta rate for the gating variable.
    """

    m_inf = a/(a + b)
    return m_inf


def calc_ah(u):
    """
    Calculates the alpha rate for the gating variable h.

    Parameters
    ----------
    u : float
        Membrane potential.
    """

    ah = calc_where(u >= -40, 0, 0.135 * exp(-(80 + u)/6.8))
    return ah


def calc_bh(u):
    """
    Calculates the beta rate for the gating variable h.

    Parameters
    ----------
    u : float
        Membrane potential.
    """

    bh = calc_where(u >= -40,
                    1 / (0.13 * (1 + exp(-(u + 10.66) / 11.1))),
                    3.56 * exp(0.079 * u) + 310000 * exp(0.35 * u)) 
    return bh


def calc_aj(u):
    """
    Calculates the alpha rate for the gating variable j.

    Parameters
    ----------
    u : float
        Membrane potential.
    """

    aj = calc_where(u >= -40, 0, 
                    ((-127140 * exp(0.2444 * u) - 0.00003474 * exp(-0.04391 * u)) * 
                     (u + 37.78)/(1 + exp(0.311*(u + 79.23)))))
    return aj


def calc_bj(u):
    """
    Calculates the beta rate for the gating variable j.

    Parameters
    ----------
    u : float
        Membrane potential.
    """

    bj = calc_where(u >= -40,
                    0.3 * exp(-0.0000002535 * u)/(1 + exp(-0.1 * (u + 32))),
                    0.1212 * exp(-0.01052 * u)/(1 + exp(-0.1378 * (u + 40.14))))
    return bj


def calc_ik1(u, gk1, ek, Cm):
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
    """

    ik1 = Cm *  gk1*(u - ek)/(1 + exp(0.07*(u + 80)))
    return ik1


def calc_tau_oa(u, kq10):
    """
    Calculates the time constant for the gating variable oa.

    Parameters
    ----------
    u : float
        Membrane potential.
    kq10 : float
        Temperature adjustment factor.
    """

    ao = 0.65/(exp(-(u + 10)/8.5) + exp(-(u - 30)/59.0))
    bo = 0.65/(2.5 + exp((u + 82)/17.0))

    tau_oa = 1/(kq10*(ao + bo))
    return tau_oa


def calc_oa_inf(u):
    """
    Calculates the steady-state value for the gating variable oa.

    Parameters
    ----------
    u : float
        Membrane potential.
    """

    oa_inf = 1/(1 + exp(-(u + 20.47)/17.54))
    return oa_inf


def calc_tau_oi(u, kq10):
    """
    Calculates the time constant for the gating variable oi.

    Parameters
    ----------
    u : float
        Membrane potential.
    kq10 : float
        Temperature adjustment factor.
    """

    aoi = 1/(18.53 + exp((u + 113.7)/10.95))
    boi = 1/(35.56 + exp(-(u + 1.26)/7.44))

    tau_oi = 1/(kq10*(aoi + boi))
    return tau_oi


def calc_oi_inf(u):
    """
    Calculates the steady-state value for the gating variable oi.

    Parameters
    ----------
    u : float
        Membrane potential.
    """

    oi_inf = 1/(1 + exp((u + 43.1)/5.3))
    return oi_inf


def calc_ito(u, oa, oi, gto, ek, Cm):
    """
    Calculates the transient outward potassium current.

    Parameters
    ----------
    u : float
        Membrane potential.
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
    """

    ito = Cm * gto*(oa**3)*oi*(u - ek)  

    return ito


def calc_tau_ua(u, kq10):
    """
    Calculates the time constant for the gating variable ua.

    Parameters
    ----------
    u : float
        Membrane potential.
    kq10 : float
        Temperature adjustment factor.
    """

    aua = 0.65/(exp(-(u + 10)/8.5) + exp(-(u - 30)/59.0))
    bua = 0.65/(2.5 + exp((u + 82)/17.0))

    tau_ua = 1/(kq10*(aua + bua))
    return tau_ua


def calc_ua_inf(u):
    """
    Calculates the steady-state value for the gating variable ua.

    Parameters
    ----------
    u : float
        Membrane potential.
    """

    ua_inf = 1/(1 + exp(-(u + 30.3)/9.6))
    return ua_inf


def calc_tau_ui(u, kq10):
    """
    Calculates the time constant for the gating variable ui.

    Parameters
    ----------
    u : float
        Membrane potential.
    kq10 : float
        Temperature adjustment factor.
    """

    aui = 1/(21 + exp(-(u - 185)/28.0))
    bui = exp((u - 158)/16.0)

    tau_ui = 1/(kq10*(aui + bui))
    return tau_ui


def calc_ui_inf(u):
    """
    Calculates the steady-state value for the gating variable ui.

    Parameters
    ----------
    u : float
        Membrane potential.
    """

    ui_inf = 1/(1 + exp((u - 99.45)/27.48))
    return ui_inf


def calc_ikur(u, ua, ui, ek, Cm):
    """
    Calculates the ultra-rapid delayed rectifier potassium current.

    Parameters
    ----------
    u : float
        Membrane potential.
    ua : float
        Gating variable ua.
    ui : float
        Gating variable ui.
    ek : float
        Equilibrium potential for potassium.
    Cm : float
        Cell membrane capacitance.
    """

    gkur = 0.005 + 0.05/(1 + exp(-(u - 15)/13.0))

    ikur = Cm * gkur*(ua**3)*ui*(u - ek)

    return ikur

def calc_tau_xr(u):
    """
    Calculates the time constant for the gating variable xr.

    Parameters
    ----------
    u : float
        Membrane potential.
    
    Note
    ----
        The singularities at u = -14.1 and u = 3.3328 are handled using
        lim x->0 : x / (1 - exp(x)) = -1
    """

    axr = calc_where(abs(u + 14.1) < 1e-5,
                     0.0003 * 5,
                     0.0003 * (u + 14.1) / (1 - exp(-(u + 14.1) / 5)))
    bxr = calc_where(abs(u - 3.3328) < 1e-5,
                     0.000073898 * 5.1237,
                     0.000073898 * (u - 3.3328) / (exp((u - 3.3328) / 5.1237) - 1))

    tau_xr = 1/(axr + bxr)
    return tau_xr


def calc_xr_inf(u):
    """
    Calculates the steady-state value for the gating variable xr.

    Parameters
    ----------
    u : float
        Membrane potential.
    """

    xr_inf = 1/(1 + exp(-(u + 14.1)/6.5))
    return xr_inf


def calc_ikr(u, xr, gkr, ek, Cm):
    """
    Calculates the rapid delayed rectifier potassium current.

    Parameters
    ----------
    u : float
        Membrane potential.
    xr : float
        Gating variable xr.
    gkr : float
        Maximum conductance for the rapid delayed rectifier potassium current.
    ek : float
        Equilibrium potential for potassium.
    Cm : float
        Cell membrane capacitance.
    """
    gkr = 0.0294 # * sqrt(ko / 5.4)

    ikr = Cm * (gkr*xr*(u - ek))/(1 + exp((u + 15)/22.4))

    return ikr


def calc_tau_xs(u):
    """
    Calculates the time constant for the gating variable xs.

    Parameters
    ----------
    u : float
        Membrane potential.
    
    Note
    ----
        The singularity at u = 19.9 mV is handled using
        lim x->0 : x / (1 - exp(x)) = -1
    """

    axs = calc_where(abs(u - 19.9) < 1e-5,
                     0.00004 * 17,
                     0.00004 * (u - 19.9) / (1 - exp(-(u - 19.9)/17)))
    bxs = calc_where(abs(u - 19.9) < 1e-5, 
                     0.000035 * 9,
                     0.000035 * (u - 19.9) / (exp((u - 19.9)/9) - 1))

    tau_xs = 1/(2*(axs + bxs))
    return tau_xs


def calc_xs_inf(u):
    """
    Calculates the steady-state value for the gating variable xs.

    Parameters
    ----------
    u : float
        Membrane potential.
    """

    xs_inf = 1 / sqrt(1 + exp(-(u - 19.9)/12.7))
    return xs_inf


def calc_iks(u, xs, gks, ek, Cm):
    """
    Calculates the slow delayed rectifier potassium current.
    """
    
    iks = Cm * gks*(xs**2)*(u - ek)

    return iks


def calc_tau_d(u):
    """
    Calculates the time constant for the gating variable d.

    Parameters
    ----------
    u : float
        Membrane potential.

    Note
    ----
        The singularity at u = -10 mV is handled using Taylor expansion of
        (exp(x) - 1) / (exp(x) + 1) = tanh(x/2) around x = 0.
    """

    tau_d = calc_where(abs(u + 10) < 1e-5,
                       1 / (6.24 * 0.035) * (1. / 2. - 1. / 24. * ((u + 10) / 6.24) ** 2),
                       (1 - exp(-(u + 10) / 6.24)) / (0.035 *(u + 10)*(1 + exp(-(u + 10)/6.24))))
    return tau_d


def calc_d_inf(u):
    """
    Calculates the steady-state value for the gating variable d.

    Parameters
    ----------
    u : float
        Membrane potential.
    """

    d_inf = 1/(1 + exp(-(u + 10)/8.0))
    return d_inf


def calc_tau_f(u):
    """
    Calculates the time constant for the gating variable f.

    Parameters
    ----------
    u : float
        Membrane potential.
    """

    tau_f = 9/(0.0197*exp(-(0.0337**2)*((u + 10)**2)) + 0.02)
    return tau_f


def calc_f_inf(u):
    """
    Calculates the steady-state value for the gating variable f.

    Parameters
    ----------
    u : float
        Membrane potential.
    """

    f_inf = 1/(1 + exp((u + 28)/6.9))
    return f_inf


def calc_tau_fca():
    """
    Calculates the time constant for the gating variable fca.

    Parameters
    ----------
    u : float
        Membrane potential.
    """

    tau_fca = 2
    return tau_fca


def calc_fca_inf(cai):
    """
    Calculates the steady-state value for the gating variable fca.

    Parameters
    ----------
    cai : float
        Intracellular calcium concentration.
    """

    fca_inf = 1/(1 + cai/0.00035)
    return fca_inf


def calc_ical(u, d, f, gcal, fca, Cm):
    """
    Calculates the L-type calcium current.

    Parameters
    ----------
    u : float
        Membrane potential.
    d : float
        Gating variable d.
    f : float
        Gating variable f.
    gcal : float
        Maximum conductance for the L-type calcium current.
    fca : float
        Gating variable fca.
    Cm : float
        Cell membrane capacitance.      
    """

    ical = Cm * gcal*d*f*fca*(u - 65) 

    return ical


def calc_inak(inakmax, nai, nao, ko, kmnai, kmko, F, u, R, T, Cm):
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
    """

    s = (1/7.0)*(exp(nao/67.3) - 1)
    fnak = 1/(1 + 0.1245*exp(-0.1*(F*u)/(R*T)) + 0.0365*s*exp(-(F*u)/(R*T)))

    inak = Cm *  inakmax*fnak*(1/(1 + (kmnai/nai)**1.5))*(ko/(ko + kmko))

    return inak


def calc_inaca(inacamax, nai, nao, cai, cao, kmnancx, kmcancx, ksatncx, F, u, R, T, Cm):
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
    inaca = Cm * numerator / denominator
    # inaca = numerator / denominator

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


def calc_Fn(irel, ical, inaca, F, Vrel):
    """
    Calculates the sarcoplasmic Ca2+ flux signal for irel.

    Parameters
    ----------
    irel : float
        Calcium release from the JSR.
    ical : float
        L-type calcium current.
    inaca : float
        Sodium-calcium exchanger current.
    F : float
        Faraday's constant.
    Vrel : float
        Volume of the release compartment.
    """

    Fn = 1e-12 * Vrel * irel - ((5 * 1e-13) / F) * (0.5 * ical - 0.2 * inaca)
    return Fn


def calc_tau_urel():
    """
    Calculates the time constant for the gating variable urel.
    """
    tau_urel = 8
    return tau_urel


def calc_urel_inf(Fn):
    """
    Calculates the steady-state value for the gating variable urel.

    Parameters
    ----------
    Fn : float
        Function Fn used in the gating variable equations.
    irel : float
        Calcium release from the JSR.
    ical : float
        L-type calcium current.
    inaca : float
        Sodium-calcium exchanger current.
    F : float
        Faraday's constant.
    Vrel : float
        Volume of the release compartment.
    """

    urel_inf = 1/(1 + exp(-(Fn - 3.4175e-13)/13.67e-16))
    return urel_inf


def calc_tau_vrel(Fn):
    """
    Calculates the time constant for the gating variable vrel.

    Parameters
    ----------
    Fn : float
        Function Fn used in the gating variable equations.
    """

    tau_vrel = 1.91 + 2.09/(1 + exp(-(Fn - 3.4175e-13)/13.67e-16))
    return tau_vrel


def calc_vrel_inf(Fn):
    """
    Calculates the steady-state value for the gating variable vrel.

    Parameters
    ----------
    Fn : float
        Function Fn used in the gating variable equations.
    """

    vrel_inf = 1 - 1/(1 + exp(-(Fn - 6.835e-14)/13.67e-16))
    return vrel_inf


def calc_tau_wrel(u):
    """
    Calculates the time constant for the gating variable wrel.

    Parameters
    ----------
    u : float
        Membrane potential.
    """

    tau_wrel = 6 * (1 - exp(-(u - 7.9) / 5.0)) / ((1 + 0.3 * exp(-(u - 7.9) / 5.0)) * (u - 7.9))
    return tau_wrel


def calc_wrel_inf(u):
    """
    Calculates the steady-state value for the gating variable wrel.

    Parameters
    ----------
    u : float
        Membrane potential.
    """

    wrel_inf = 1 - 1/(1 + exp(-(u - 40)/17.0))
    return wrel_inf


def calc_irel(urel, vrel, wrel, krel, carel, cai): 
    """
    Calculates the calcium release from the JSR.

    Parameters
    ----------
    urel : float    
        Gating variable urel.
    vrel : float
        Gating variable vrel.
    irel : float
        Calcium release from the JSR.
    wrel : float
        Gating variable wrel.
    krel : float
        Scaling factor for the calcium release.
    carel : float
        Calcium concentration in the release compartment.
    cai : float
        Intracellular calcium concentration.
    """

    irel = krel * (urel**2) * vrel * wrel * (carel - cai)

    return irel 


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