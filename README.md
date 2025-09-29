## Courtemanche finitewave model

This model describes the ionic currents and action potential dynamics of human atrial myocytes. It includes detailed formulations for major ionic currents (fast sodium current, L-type calcium current, inward rectifier potassium current, transient outward potassium current, rapid and slow delayed rectifier potassium currents, and Na⁺/Ca²⁺ exchanger), as well as calcium handling mechanisms.

The Courtemanche model is widely used as a reference atrial electrophysiology model. It has served as the basis for many subsequent atrial modeling studies, including investigations of atrial fibrillation and drug effects.

### Reference
Courtemanche M, Ramirez RJ, Nattel S. Ionic mechanisms underlying human atrial action potential properties: insights from a mathematical model. Am J Physiol. 1998 Jul;275(1):H301-21. 

DOI: 10.1152/ajpheart.1998.275.1.H301

### How to use (quickstart)
```bash
python -m examples.courtemanche_example
```

### How to test
```bash
python -m pytest -q
```

### Repository structure
```text
.
├── courtemanche/                    # equations package (ops.py)
│   ├── __init__.py
│   └── ops.py                       # model equations (pure functions)
├── implementation/                  # 0D model implementation
│   ├── __init__.py
│   └── courtemanche_0d.py
├── example/
│   └── courtemanche_example.py      # minimal script to run a short trace
├── tests/
│   └── test.py                      # smoke test; reproducibility checks
├── .gitignore
├── LICENSE                          # MIT
├── pyproject.toml                   # configuration file
└── README.md                        # this file
```

### Variables
- `u = -84.5` - Membrane potential (mV)
- `nai = 11.2` - intracellular Na⁺ concentration (mM)
- `ki = 139` - intracellular K⁺ concentration (mM)
- `cai = 0.000102` - intracellular Ca²⁺ concentration (mM)
- `caup = 1.6` - Ca²⁺ concentration in uptake compartment (mM)
- `carel = 1.1` - Ca²⁺ concentration in release compartment (mM)
- `m = 0.00291 ` - activation gate of fast Na⁺ current
- `h = 0.965` - fast inactivation gate of I_Na
- `j = 0.978` - slow inactivation gate of I_Na
- `d = 0.000137` - activation gate of L-type Ca²⁺ current
- `f = 0.999837` - inactivation gate of I_CaL
- `oa = 0.000592` - activation gate of transient outward K⁺ current
- `oi = 0.9992` - inactivation gate of I_to
- `ua = 0.003519` - activation gate of ultra-rapid K⁺ current
- `ui = 0.9987` - inactivation gate of I_Kur
- `xs = 0.0187` - activation gate of slow delayed rectifier K⁺ current
- `xr = 0.0000329` - activation gate of rapid delayed rectifier K⁺ current
- `fca = 0.775` - Ca²⁺-dependent inactivation gate of I_CaL
- `irel  = 0` - Ca²⁺ release current from SR
- `vrel  = 1` - SR release flux variable (voltage-dependent)
- `urel  = 0` - SR release gate variable
- `wrel  = 0.9` - SR recovery variable


### Parameters
- `gna = 7.8` - maximum conductance of fast sodium current (mS/cm²)
- `gnab = 0.000674` - maximum conductance of background sodium current (mS/cm²)
- `gk1 = 0.09` - maximum conductance of inward rectifier potassium current (mS/cm²)
- `gkr = 0.0294` - maximum conductance of rapid delayed rectifier potassium current (mS/cm²)
- `gks = 0.129` - maximum conductance of slow delayed rectifier potassium current (mS/cm²)
- `gto = 0.1652` - maximum conductance of transient outward potassium current (mS/cm²)
- `gcal = 0.1238` - maximum conductance of L-type calcium current (mS/cm²)
- `gcab = 0.00113` - maximum conductance of background calcium current (mS/cm²)

- `Vc = 20100` - cell volume (µL)
- `Vj = Vc * 0.68` - junction volume (µL)
- `Vup = Vj * 0.06 * 0.92` - uptake compartment volume (µL)
- `Vrel = Vj * 0.06 * 0.08` - release compartment volume (µL)

- `ibk = 0.0` - background current (µA/cm²)
- `cao = 1.8` - extracellular Ca²⁺ concentration (mM)
- `nao = 140` - extracellular Na⁺ concentration (mM)
- `ko = 5.4` - extracellular K⁺ concentration (mM)

- `caupmax = 15` - maximum Ca²⁺ uptake rate (mM/ms)
- `kup = 0.00092` - uptake rate constant (mM)

- `kmnai = 10` - Michaelis-Menten constant for Na⁺ (mM)
- `kmko = 1.5` - Michaelis-Menten constant for K⁺ (mM)
- `kmnancx = 87.5` - Michaelis-Menten constant for Na⁺/Ca²⁺ exchanger (mM)
- `kmcancx = 1.38` - Michaelis-Menten constant for Ca²⁺/Na⁺ exchanger (mM)
- `ksatncx = 0.1` - saturation constant for Na⁺/Ca²⁺ exchanger

- `kmcmdn = 0.00238` - Michaelis-Menten constant for calmodulin (mM)
- `kmtrpn = 0.0005` - Michaelis-Menten constant for troponin (mM)
- `kmcsqn = 0.8` - Michaelis-Menten constant for calsequestrin (mM)

- `trpnmax = 0.07` - maximum troponin concentration (mM)
- `cmdnmax = 0.05` - maximum calmodulin concentration (mM)
- `csqnmax = 10.0` - maximum calsequestrin concentration (mM)
- `inacamax = 1600` - maximum Na⁺ current (µA/cm²)
- `inakmax = 0.6` - maximum Na⁺/K⁺ pump current (µA/cm²)
- `ipcamax = 0.275` - maximum Ca²⁺ pump current (µA/cm²)
- `krel = 30` - release rate constant

- `iupmax = 0.005` - maximum uptake current (µA/cm²)

- `kq10 = 3` - temperature sensitivity coefficient

