#
# Compare lithium-ion battery models
#
import pybamm

pybamm.set_logging_level("INFO")

# load models
models = [
    pybamm.lithium_ion.SPM({"loss of active material": "stress-driven"}),
    pybamm.lithium_ion.SPMe({"loss of active material": "stress-driven"}),
    pybamm.lithium_ion.DFN({"loss of active material": "stress-driven"}),
    pybamm.lithium_ion.NewmanTobias({"loss of active material": "stress-driven"}),
]

# create and run simulations
sims = []
for model in models:
    sim = pybamm.Simulation(model)
    sim.solve([0, 3600])
    sims.append(sim)

# plot
pybamm.dynamic_plot(sims)
