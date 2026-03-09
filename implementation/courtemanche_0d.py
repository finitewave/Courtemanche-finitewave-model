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
        ena, ek, eca = ops.calc_equilibrum_potentials(self.variables["nai"], self.parameters["nao"], self.variables["ki"],
                                                      self.parameters["ko"], self.variables["cai"], self.parameters["cao"], self.parameters["R"], self.parameters["T"], self.parameters["F"])
        
        var_m = ops.calc_gating_m(self.variables["m"], self.variables["u"], self.dt)
        var_h = ops.calc_gating_h(self.variables["h"], self.variables["u"], self.dt)
        var_j = ops.calc_gating_j(self.variables["j"], self.variables["u"], self.dt)

        ina = ops.calc_ina(self.variables["u"], self.variables["m"], self.variables["h"], self.variables["j"], self.parameters["gna"], ena, self.parameters["Cm"])
        ik1 = ops.calc_ik1(self.variables["u"], self.parameters["gk1"], ek, self.parameters["Cm"])
        ito, var_oa, var_oi = ops.calc_ito(self.variables["u"], self.dt, self.parameters["kq10"], self.variables["oa"], 
                                                                       self.variables["oi"], self.parameters["gto"], ek, self.parameters["Cm"])
        ikur, var_ua, var_ui = ops.calc_ikur(self.variables["u"], self.dt, self.parameters["kq10"], self.variables["ua"], 
                                                                         self.variables["ui"], ek, self.parameters["Cm"])
        ikr, var_xr = ops.calc_ikr(self.variables["u"], self.dt, self.variables["xr"], self.parameters["gkr"], ek, self.parameters["Cm"])
        iks, var_xs = ops.calc_iks(self.variables["u"], self.dt, self.variables["xs"], self.parameters["gks"], ek, self.parameters["Cm"])
        ical, var_d, var_f, var_fca = ops.calc_ical(self.variables["u"], self.dt, self.variables["d"], self.variables["f"], 
                                                                                              self.variables["cai"], self.parameters["gcal"], self.variables["fca"], self.parameters["Cm"])
        inak = ops.calc_inak(self.parameters["inakmax"], self.variables["nai"], self.parameters["nao"], self.parameters["ko"], self.parameters["kmnai"], 
                             self.parameters["kmko"], self.parameters["F"], self.variables["u"], self.parameters["R"], self.parameters["T"], self.parameters["Cm"])
        inaca = ops.calc_inaca(self.parameters["inacamax"], self.variables["nai"], self.parameters["nao"], self.variables["cai"], self.parameters["cao"], 
                               self.parameters["kmnancx"], self.parameters["kmcancx"], self.parameters["ksatncx"], self.parameters["F"], self.variables["u"], 
                               self.parameters["R"], self.parameters["T"], self.parameters["Cm"])
        ibca = ops.calc_ibca(self.parameters["gcab"], eca, self.variables["u"], self.parameters["Cm"])
        ibna = ops.calc_ibna(self.parameters["gnab"], ena, self.variables["u"], self.parameters["Cm"])
        ipca = ops.calc_ipca(self.parameters["ipcamax"], self.variables["cai"], self.parameters["Cm"])
        var_irel, var_urel, var_vrel, var_wrel = ops.calc_irel(self.dt, self.variables["urel"], self.variables["vrel"], self.variables["irel"], self.variables["wrel"], ical, inaca, self.parameters["krel"], self.variables["carel"], self.variables["cai"], self.variables["u"], self.parameters["F"], self.parameters["Vrel"])
        itr = ops.calc_itr(self.variables["caup"], self.variables["carel"])
        iup = ops.calc_iup(self.parameters["iupmax"], self.variables["cai"], self.parameters["kup"])
        iupleak = ops.calc_iupleak(self.variables["caup"], self.parameters["caupmax"], self.parameters["iupmax"])

        var_caup = ops.calc_dcaup(iup, iupleak, itr, self.parameters["Vrel"], self.parameters["Vup"])
        var_nai = ops.calc_dnai(inak, inaca, ibna, ina, self.parameters["F"], self.parameters["Vj"])

        var_ki = ops.calc_dki(inak, ik1, ito, ikur, ikr, iks, self.parameters["ibk"], self.parameters["F"], self.parameters["Vj"])
        var_cai = ops.calc_dcai(self.variables["cai"], inaca, ipca, ical, ibca, iup, iupleak, self.variables["irel"], self.parameters["Vrel"], 
                                                      self.parameters["Vup"], self.parameters["trpnmax"], self.parameters["kmtrpn"], self.parameters["cmdnmax"], self.parameters["kmcmdn"], 
                                                      self.parameters["F"], self.parameters["Vj"])

        var_carel = ops.calc_dcarel(self.variables["carel"], itr, self.variables["irel"], self.parameters["csqnmax"], self.parameters["kmcsqn"])

        # Update variables:
        self.variables["m"] = var_m
        self.variables["h"] = var_h
        self.variables["j"] = var_j
        self.variables["oa"] = var_oa
        self.variables["oi"] = var_oi
        self.variables["ua"] = var_ua
        self.variables["ui"] = var_ui
        self.variables["xr"] = var_xr
        self.variables["xs"] = var_xs
        self.variables["d"] = var_d
        self.variables["f"] = var_f
        self.variables["fca"] = var_fca
        self.variables["urel"] = var_urel
        self.variables["vrel"] = var_vrel
        self.variables["irel"] = var_irel
        self.variables["wrel"] = var_wrel
        self.variables["caup"] += self.dt*var_caup
        self.variables["nai"] += self.dt*var_nai
        self.variables["ki"] += self.dt*var_ki
        self.variables["cai"] += self.dt*var_cai
        self.variables["carel"] += self.dt*var_carel

        self.variables["u"] += self.dt*(-ops.calc_rhs(ina, ik1, ito, ikur, ikr, iks, ical, ipca, inak, inaca, ibna, ibca, self.parameters["Cm"]) + sum(stim.stim(t=self.dt*i) for stim in self.stimulations))

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