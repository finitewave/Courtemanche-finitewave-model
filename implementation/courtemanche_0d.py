"""
This module provides a simple interface to run the model in a 0D setting,
i.e., without spatial dimensions. It includes class for defining stimulation protocols
and a class for the 0D model itself.

"""

from courtemanche import ops


class Stimulation:
    """
    Stimulus protocol for the 0D model.

    Parameters
    ----------
    t_start : float
        Start time (ms) of the first stimulus window.
    duration : float
        Duration (ms) of a single pulse.
    amplitude : float
        Pulse amplitude in the same units as du/dt contribution (typically "units/ms").

    Method
    ------
    stim(t: float) -> float
        Returns the instantaneous stimulus value at time t.

    """

    def __init__(self, t_start: float, duration: float, amplitude: float):
        self.t_start = t_start
        self.duration = duration
        self.amplitude = amplitude

    def stim(self, t: float) -> float:
        return self.amplitude if self.t_start <= t < self.t_start + self.duration else 0.0


class Courtemanche0D:
    """
    Courtemanche OD implementation.

    Parameters
    ----------

    dt : float
        Time step size (ms).
    stimulations : list[Stimulation]
        List of stimulation protocols to apply during the simulation.

    Attributes
    ----------
    variables : dict[str, float]
        Current state variables of the model.
    parameters : dict[str, float]
        Model parameters.
    history : dict[str, list[float]]
        Time history of state variables for post-processing.
    
    Methods
    -------
    step(i: int)
        Perform a single time step update.
    run(t_max: float)
        Run the simulation up to time t_max.
    """
    def __init__(self, dt: float, stimulations: list[Stimulation]):
        self.dt = dt
        self.stimulations = stimulations
        self.variables = ops.get_variables()
        self.parameters = ops.get_parameters()
        self.history = {s: [] for s in self.variables}

    def step(self, i: int):
        """
        Perform a single time step update.

        Parameters
        ----------
        i : int
            Current time step index.
        """
        u_old = self.variables["u"]
        m_old = self.variables["m"]
        h_old = self.variables["h"]
        j_old = self.variables["j"]

        oa_old = self.variables["oa"]
        oi_old = self.variables["oi"]
        ua_old = self.variables["ua"]
        ui_old = self.variables["ui"]
        xr_old = self.variables["xr"]
        xs_old = self.variables["xs"]

        d_old = self.variables["d"]
        f_old = self.variables["f"]
        fca_old = self.variables["fca"]

        urel_old = self.variables["urel"]
        vrel_old = self.variables["vrel"]
        irel_old = self.variables["irel"]
        wrel_old = self.variables["wrel"]

        caup_old = self.variables["caup"]
        carel_old = self.variables["carel"]
        cai_old = self.variables["cai"]
        nai_old = self.variables["nai"]
        ki_old = self.variables["ki"]


        ena = ops.calc_ena(nai_old, self.parameters["nao"], self.parameters["R"], self.parameters["T"], self.parameters["F"])
        ek = ops.calc_ek(ki_old, self.parameters["ko"], self.parameters["R"], self.parameters["T"], self.parameters["F"])
        eca = ops.calc_eca(cai_old, self.parameters["cao"], self.parameters["R"], self.parameters["T"], self.parameters["F"])

        am = ops.calc_am(u_old)
        bm = ops.calc_bm(u_old)
        tau_m = ops.calc_tau(am, bm)
        m_inf = ops.calc_inf(am, bm)
        m_new = ops.calc_gating_variable_rush_larsen(m_old, m_inf, tau_m, self.dt)

        ah = ops.calc_ah(u_old)
        bh = ops.calc_bh(u_old) 
        tau_h = ops.calc_tau(ah, bh)
        h_inf = ops.calc_inf(ah, bh)
        h_new = ops.calc_gating_variable_rush_larsen(h_old, h_inf, tau_h, self.dt)

        aj = ops.calc_aj(u_old)
        bj = ops.calc_bj(u_old)
        tau_j = ops.calc_tau(aj, bj)
        j_inf = ops.calc_inf(aj, bj)
        j_new = ops.calc_gating_variable_rush_larsen(j_old, j_inf, tau_j, self.dt)

        ina = ops.calc_ina(
            u_old, m_old, h_old, j_old,
            self.parameters["gna"], ena, self.parameters["Cm"]
        )

        ik1 = ops.calc_ik1(
            u_old, self.parameters["gk1"], ek, self.parameters["Cm"]
        )

        tau_oa = ops.calc_tau_oa(u_old, self.parameters["kq10"])
        oa_inf = ops.calc_oa_inf(u_old)
        oa_new = ops.calc_gating_variable_rush_larsen(oa_old, oa_inf, tau_oa, self.dt)

        tau_oi = ops.calc_tau_oi(u_old, self.parameters["kq10"])
        oi_inf = ops.calc_oi_inf(u_old)
        oi_new = ops.calc_gating_variable_rush_larsen(oi_old, oi_inf, tau_oi, self.dt)

        ito = ops.calc_ito(
            u_old, oa_old, oi_old,
            self.parameters["gto"], ek, self.parameters["Cm"]
        )

        tau_ua = ops.calc_tau_ua(u_old, self.parameters["kq10"])
        ua_inf = ops.calc_ua_inf(u_old)
        ua_new = ops.calc_gating_variable_rush_larsen(ua_old, ua_inf, tau_ua, self.dt)

        tau_ui = ops.calc_tau_ui(u_old, self.parameters["kq10"])
        ui_inf = ops.calc_ui_inf(u_old)
        ui_new = ops.calc_gating_variable_rush_larsen(ui_old, ui_inf, tau_ui, self.dt)

        ikur = ops.calc_ikur(
            u_old, ua_old, ui_old,
            ek, self.parameters["Cm"]
        )

        tau_xr = ops.calc_tau_xr(u_old)
        xr_inf = ops.calc_xr_inf(u_old)
        xr_new = ops.calc_gating_variable_rush_larsen(xr_old, xr_inf, tau_xr, self.dt)

        ikr = ops.calc_ikr(
            u_old, xr_old,
            self.parameters["gkr"], ek, self.parameters["Cm"]
        )

        tau_xs = ops.calc_tau_xs(u_old)
        xs_inf = ops.calc_xs_inf(u_old)
        xs_new = ops.calc_gating_variable_rush_larsen(xs_old, xs_inf, tau_xs, self.dt)

        iks = ops.calc_iks(
            u_old, xs_old,
            self.parameters["gks"], ek, self.parameters["Cm"]
        )

        tau_d = ops.calc_tau_d(u_old)
        d_inf = ops.calc_d_inf(u_old)
        d_new = ops.calc_gating_variable_rush_larsen(d_old, d_inf, tau_d, self.dt)

        tau_f = ops.calc_tau_f(u_old)
        f_inf = ops.calc_f_inf(u_old)
        f_new = ops.calc_gating_variable_rush_larsen(f_old, f_inf, tau_f, self.dt)

        tau_fca = ops.calc_tau_fca()
        fca_inf = ops.calc_fca_inf(cai_old)
        fca_new = ops.calc_gating_variable_rush_larsen(fca_old, fca_inf, tau_fca, self.dt)

        ical = ops.calc_ical(
            u_old, d_old, f_old,
            self.parameters["gcal"], fca_old, self.parameters["Cm"]
        )

        inak = ops.calc_inak(
            self.parameters["inakmax"],
            nai_old, self.parameters["nao"],
            self.parameters["ko"], self.parameters["kmnai"],
            self.parameters["kmko"], self.parameters["F"],
            u_old, self.parameters["R"], self.parameters["T"],
            self.parameters["Cm"]
        )

        inaca = ops.calc_inaca(
            self.parameters["inacamax"],
            nai_old, self.parameters["nao"],
            cai_old, self.parameters["cao"],
            self.parameters["kmnancx"], self.parameters["kmcancx"],
            self.parameters["ksatncx"], self.parameters["F"],
            u_old, self.parameters["R"], self.parameters["T"],
            self.parameters["Cm"]
        )

        ibca = ops.calc_ibca(
            self.parameters["gcab"], eca, u_old, self.parameters["Cm"]
        )

        ibna = ops.calc_ibna(
            self.parameters["gnab"], ena, u_old, self.parameters["Cm"]
        )

        ipca = ops.calc_ipca(
            self.parameters["ipcamax"], cai_old, self.parameters["Cm"]
        )

        Fn = ops.calc_Fn(irel_old, ical, inaca, self.parameters["F"], self.parameters["Vrel"])

        tau_urel = ops.calc_tau_urel()
        urel_inf = ops.calc_urel_inf(Fn)
        urel_new = ops.calc_gating_variable_rush_larsen(urel_old, urel_inf, tau_urel, self.dt)

        tau_vrel = ops.calc_tau_vrel(Fn)
        vrel_inf = ops.calc_vrel_inf(Fn)
        vrel_new = ops.calc_gating_variable_rush_larsen(vrel_old, vrel_inf, tau_vrel, self.dt)

        tau_wrel = ops.calc_tau_wrel(u_old)
        wrel_inf = ops.calc_wrel_inf(u_old)
        wrel_new = ops.calc_gating_variable_rush_larsen(wrel_old, wrel_inf, tau_wrel, self.dt)

        irel_new = ops.calc_irel(
            urel_old, vrel_old, irel_old, wrel_old,
            self.parameters["krel"],
            carel_old, cai_old
        )

        itr = ops.calc_itr(caup_old, carel_old)
        iup = ops.calc_iup(self.parameters["iupmax"], cai_old, self.parameters["kup"])
        iupleak = ops.calc_iupleak(caup_old, self.parameters["caupmax"], self.parameters["iupmax"])

        dcaup = ops.calc_dcaup(
            iup, iupleak, itr,
            self.parameters["Vrel"], self.parameters["Vup"]
        )

        dnai = ops.calc_dnai(
            inak, inaca, ibna, ina,
            self.parameters["F"], self.parameters["Vj"]
        )

        dki = ops.calc_dki(
            inak, ik1, ito, ikur, ikr, iks,
            self.parameters["ibk"], self.parameters["F"], self.parameters["Vj"]
        )

        dcai = ops.calc_dcai(
            cai_old, inaca, ipca, ical, ibca,
            iup, iupleak, irel_old,
            self.parameters["Vrel"], self.parameters["Vup"],
            self.parameters["trpnmax"], self.parameters["kmtrpn"],
            self.parameters["cmdnmax"], self.parameters["kmcmdn"],
            self.parameters["F"], self.parameters["Vj"]
        )

        dcarel = ops.calc_dcarel(
            carel_old, itr, irel_old,
            self.parameters["csqnmax"], self.parameters["kmcsqn"]
        )

        stim_current = sum(stim.stim(t=self.dt * i) for stim in self.stimulations)

        du = -ops.calc_rhs(
            ina, ik1, ito, ikur, ikr, iks,
            ical, ipca, inak, inaca, ibna, ibca,
            self.parameters["Cm"]
        ) + stim_current

        caup_new = caup_old + self.dt * dcaup
        nai_new = nai_old + self.dt * dnai
        ki_new = ki_old + self.dt * dki
        cai_new = cai_old + self.dt * dcai
        carel_new = carel_old + self.dt * dcarel
        u_new = u_old + self.dt * du


        self.variables["m"] = m_new
        self.variables["h"] = h_new
        self.variables["j"] = j_new

        self.variables["oa"] = oa_new
        self.variables["oi"] = oi_new
        self.variables["ua"] = ua_new
        self.variables["ui"] = ui_new
        self.variables["xr"] = xr_new
        self.variables["xs"] = xs_new

        self.variables["d"] = d_new
        self.variables["f"] = f_new
        self.variables["fca"] = fca_new

        self.variables["urel"] = urel_new
        self.variables["vrel"] = vrel_new
        self.variables["irel"] = irel_new
        self.variables["wrel"] = wrel_new

        self.variables["caup"] = caup_new
        self.variables["nai"] = nai_new
        self.variables["ki"] = ki_new
        self.variables["cai"] = cai_new
        self.variables["carel"] = carel_new

        self.variables["u"] = u_new

    def run(self, t_max: float):
        """
        Run the simulation up to time t_max.
        
        Parameters
        ----------
        t_max : float
            Maximum simulation time.
        """
        n_steps = int(round(t_max/self.dt))
        for i in range(n_steps):
            self.step(i)
            for s in self.variables:
                self.history[s].append(self.variables[s])