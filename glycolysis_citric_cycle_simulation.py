import numpy as np
import matplotlib.pyplot as plt

# Initial metabolite concentrations for glycolysis(mM)
glucose = 10.0  # initial concentration of glucose
atp = 5.0  # initial concentration of ATP
adp = 2.0  # initial concentration of ADP
nadh = 0.5  # initial concentration of NADH
nad = 0.5  # initial concentration of NAD+
pyruvate = 0.0  # initial concentration of pyruvate

# Initial metabolite concentrations for the citric acid cycle (mM)
acetyl_coa = 2.0  # initial concentration of acetyl-CoA
oxaloacetate = 2.0  # initial concentration of oxaloacetate
citrate = 0.0  # initial concentration of citrate
atp_citric = 15.0  # initial concentration of ATP
adp_citric = 2.0  # initial concentration of ADP
nadh_citric = 0.5  # initial concentration of NADH
nad_citric = 0.5  # initial concentration of NAD+

# Time parameters
dt = 0.1  # time step (s)
sim_time = 50.0  # total simulation time (s)
num_steps = int(sim_time / dt)

# Arrays to store metabolite concentrations over time for glycolysis
glucose_conc = np.zeros(num_steps)
atp_conc = np.zeros(num_steps)
adp_conc = np.zeros(num_steps)
nadh_conc = np.zeros(num_steps)
nad_conc = np.zeros(num_steps)
pyruvate_conc = np.zeros(num_steps)

# Arrays to store metabolite concentrations over time for the citric acid cycle
acetyl_coa_conc = np.zeros(num_steps)
oxaloacetate_conc = np.zeros(num_steps)
citrate_conc = np.zeros(num_steps)
atp_citric_conc = np.zeros(num_steps)
adp_citric_conc = np.zeros(num_steps)
nadh_citric_conc = np.zeros(num_steps)
nad_citric_conc = np.zeros(num_steps)


# Define rate equations for each step in glycolysis
def glycolysis_reaction1(glucose_, atp):
    vmax = 1.0
    km_glucose = 1.0
    km_atp = 0.5
    return vmax * glucose_ / (km_glucose + glucose_) * atp / (km_atp + atp)


def glycolysis_reaction2(glucose):
    vmax = 1.0
    km_glucose = 1.0
    return vmax * glucose / (km_glucose + glucose)


def glycolysis_reaction3(glucose):
    vmax = 1.0
    km_glucose = 1.0
    return vmax * glucose / (km_glucose + glucose)


# Define rate equations for each step in the citric acid cycle
def citric_acid_cycle_reaction1(acetyl_COA, oxaloacetate, atp):
    vmax = 1.0
    km_acetyl_coa = 1.0
    km_oxaloacetate = 0.5
    km_atp = 0.5
    return vmax * acetyl_COA * oxaloacetate / (
            (km_acetyl_coa + acetyl_COA) * (km_oxaloacetate + oxaloacetate) * (km_atp + atp))


def citric_acid_cycle_reaction2(citrate):
    vmax = 1.0
    km_citrate = 1.0
    return vmax * citrate / (km_citrate + citrate)


# Run simulations
for i in range(num_steps):
    # Glycolysis reactions
    reaction_rate1 = glycolysis_reaction1(glucose, atp)
    reaction_rate2 = glycolysis_reaction2(glucose)
    reaction_rate3 = glycolysis_reaction3(glucose)

    glucose -= reaction_rate1 * dt
    atp -= reaction_rate1 * dt
    adp += reaction_rate1 * dt
    nad += reaction_rate1 * dt
    nadh -= reaction_rate1 * dt
    pyruvate += reaction_rate3 * dt

    glucose_conc[i] = glucose
    atp_conc[i] = atp
    adp_conc[i] = adp
    nad_conc[i] = nad
    nadh_conc[i] = nadh
    pyruvate_conc[i] = pyruvate

    # Citric acid cycle reactions
    reaction_rate1_citric = citric_acid_cycle_reaction1(acetyl_coa, oxaloacetate, atp_citric)
    reaction_rate2_citric = citric_acid_cycle_reaction2(citrate)

    acetyl_coa -= reaction_rate1_citric * dt
    oxaloacetate -= reaction_rate1_citric * dt
    citrate += reaction_rate1_citric * dt
    atp_citric -= reaction_rate1_citric * dt
    adp_citric += reaction_rate1_citric * dt
    nad_citric += reaction_rate1_citric * dt
    nadh_citric -= reaction_rate1_citric * dt

    citrate -= reaction_rate2_citric * dt
    acetyl_coa += reaction_rate2_citric * dt
    oxaloacetate += reaction_rate2_citric * dt
    atp_citric += reaction_rate2_citric * dt
    adp_citric -= reaction_rate2_citric * dt
    nad_citric -= reaction_rate2_citric * dt
    nadh_citric += reaction_rate2_citric * dt

    acetyl_coa_conc[i] = acetyl_coa
    oxaloacetate_conc[i] = oxaloacetate
    citrate_conc[i] = citrate
    atp_citric_conc[i] = atp_citric
    adp_citric_conc[i] = adp_citric
    nad_citric_conc[i] = nad_citric
    nadh_citric_conc[i] = nadh_citric

    # Print values at specific time points
    if i == num_steps // 4:
        print("\nAt 25% of simulation time:")
        print("Glycolysis - Glucose concentration:", glucose)
        print("Glycolysis - ATP concentration:", atp)
        print("Glycolysis - Pyruvate concentration:", pyruvate)
        print("Citric Acid Cycle - Acetyl-CoA concentration:", acetyl_coa)
        print("Citric Acid Cycle - Oxaloacetate concentration:", oxaloacetate)
        print("Citric Acid Cycle - Citrate concentration:", citrate)
    elif i == num_steps // 2:
        print("\nAt 50% of simulation time:")
        print("Glycolysis - Glucose concentration:", glucose)
        print("Glycolysis - ATP concentration:", atp)
        print("Glycolysis - Pyruvate concentration:", pyruvate)
        print("Citric Acid Cycle - Acetyl-CoA concentration:", acetyl_coa)
        print("Citric Acid Cycle - Oxaloacetate concentration:", oxaloacetate)
        print("Citric Acid Cycle - Citrate concentration:", citrate)
    elif i == 3 * (num_steps // 4):
        print("\nAt 75% of simulation time:")
        print("Glycolysis - Glucose concentration:", glucose)
        print("Glycolysis - ATP concentration:", atp)
        print("Glycolysis - Pyruvate concentration:", pyruvate)
        print("Citric Acid Cycle - Acetyl-CoA concentration:", acetyl_coa)
        print("Citric Acid Cycle - Oxaloacetate concentration:", oxaloacetate)
        print("Citric Acid Cycle - Citrate concentration:", citrate)

# Plot results
time = np.linspace(0, sim_time, num_steps)

plt.plot(time, pyruvate_conc, label='Pyruvate')
plt.xlabel('Time (s)')
plt.ylabel('Concentration (mM)')
plt.title('Glycolysis')
plt.legend()

plt.subplot(2, 2, 2)
plt.plot(time, acetyl_coa_conc, label='Acetyl-CoA')
plt.plot(time, oxaloacetate_conc, label='Oxaloacetate')
plt.plot(time, citrate_conc, label='Citrate')
plt.xlabel('Time (s)')
plt.ylabel('Concentration (mM)')
plt.title('Citric Acid Cycle')
plt.legend()

plt.subplot(2, 2, 3)
plt.plot(time, atp_conc, label='ATP')
plt.plot(time, adp_conc, label='ADP')
plt.xlabel('Time (s)')
plt.ylabel('Concentration (mM)')
plt.title('ATP/ADP Concentration')
plt.legend()

plt.subplot(2, 2, 4)
plt.plot(time, nadh_conc, label='NADH')
plt.plot(time, nad_conc, label='NAD+')
plt.xlabel('Time (s)')
plt.ylabel('Concentration (mM)')
plt.title('NADH/NAD+ Concentration')
plt.legend()

plt.tight_layout()
plt.show()
