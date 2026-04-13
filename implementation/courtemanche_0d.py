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
        self.stim_history = []

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

        res = ops.ionic_step(
            self.dt, u_old, m_old, h_old, j_old,
            oa_old, oi_old, ua_old, ui_old, xr_old, xs_old,
            d_old, f_old, fca_old,
            urel_old, vrel_old, irel_old, wrel_old,
            caup_old, carel_old, cai_old, nai_old, ki_old,
            self.parameters["R"], self.parameters["T"], self.parameters["F"],
            self.parameters["Cm"], self.parameters["Vc"], self.parameters["Vj"],
            self.parameters["Vup"], self.parameters["Vrel"],
            self.parameters["ko"], self.parameters["nao"], self.parameters["cao"],
            self.parameters["gna"], self.parameters["gk1"], self.parameters["gto"],
            self.parameters["gkr"], self.parameters["gks"], self.parameters["gcal"],
            self.parameters["gcab"], self.parameters["gnab"],
            self.parameters["inakmax"], self.parameters["inacamax"], self.parameters["ipcamax"],
            self.parameters["iupmax"], self.parameters["kq10"], self.parameters["gamma"],
            self.parameters["kmnai"], self.parameters["kmko"], self.parameters["kmnancx"],
            self.parameters["kmcancx"], self.parameters["ksatncx"], self.parameters["krel"],
            self.parameters["kup"], self.parameters["caupmax"], self.parameters["cmdnmax"],
            self.parameters["trpnmax"], self.parameters["csqnmax"], self.parameters["kmcmdn"],
            self.parameters["kmtrpn"], self.parameters["kmcsqn"], self.parameters["ibk"],
        )
        (rhs, m_new, h_new, j_new, oa_new, oi_new, ua_new, ui_new, xr_new, xs_new,
         d_new, f_new, fca_new, urel_new, vrel_new, irel_new, wrel_new,
         caup_new, carel_new, cai_new, nai_new, ki_new) = res

        
        stim_curr = self.dt * sum(stim.stim(t=self.dt*i) for stim in self.stimulations)
        self.stim_history.append(stim_curr)

        self.variables["u"] += self.dt * rhs + stim_curr  

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