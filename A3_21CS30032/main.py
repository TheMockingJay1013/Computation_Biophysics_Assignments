
# get the fragrment files from the robetta fragment server
import requests

url = 'http://old.robetta.org/downloads/fragments/82163/aat000_03_05.200_v1_3'
output_file = 'aat000_03.frag'

response = requests.get(url)

if response.status_code == 200:
    with open(output_file, 'wb') as file:
        file.write(response.content)


url = 'http://old.robetta.org/downloads/fragments/82163/aat000_09_05.200_v1_3'
output_file = 'aat000_09.frag'

response = requests.get(url)

if response.status_code == 200:
    with open(output_file, 'wb') as file:
        file.write(response.content)

import pyrosetta
from pyrosetta import rosetta


# Initialization
pyrosetta.init("-in::file::fullatom -mute all")

# Create starting pose from sequence
# sequence = "MLSDEDFKAFGMTRSAFANLPLWKQQNLKKEKLLF"
sequence = "MLSDEDFKAVFGMTRSAFANLPLWKQQNLKKEKGLF"
pose = pyrosetta.pose_from_sequence(sequence, "fa_standard")


# Linearize pose (extended conformation)
for i in range(1, pose.total_residue() + 1):
    pose.set_phi(i, -150)
    pose.set_psi(i, 150)
    pose.set_omega(i, 180)

# Convert to centroid mode
to_centroid = rosetta.protocols.simple_moves.SwitchResidueTypeSetMover("centroid")
to_centroid.apply(pose)

# Setup MoveMap
movemap = rosetta.core.kinematics.MoveMap()
movemap.set_bb(True)

# Load fragment files
frag9 = rosetta.core.fragment.ConstantLengthFragSet(9,"aat000_09.frag")
frag3 = rosetta.core.fragment.ConstantLengthFragSet(3,"aat000_03.frag")


# Create fragment movers
fragmover9 = rosetta.protocols.simple_moves.ClassicFragmentMover(frag9, movemap)
fragmover3 = rosetta.protocols.simple_moves.ClassicFragmentMover(frag3, movemap)
# fragmover9.set_check_ss(False)
# fragmover3.set_check_ss(False)

# Repeat movers and sequence mover
repeated_9mer = rosetta.protocols.moves.RepeatMover(fragmover9, 3)
repeated_3mer = rosetta.protocols.moves.RepeatMover(fragmover3, 9)
seq_mover = rosetta.protocols.moves.SequenceMover()
seq_mover.add_mover(repeated_9mer)
seq_mover.add_mover(repeated_3mer)

# Setup Monte Carlo
score3 = rosetta.core.scoring.ScoreFunctionFactory.create_score_function("score3")
kT = 3.0
mc = rosetta.protocols.moves.MonteCarlo(pose, score3, kT)
trial_mover = rosetta.protocols.moves.TrialMover(seq_mover, mc)

score_list = []  # list to store scores
cycles = 50000  # number of cycles for the folding simulation


# Run the folding simulation
for i in range(cycles):
    trial_mover.apply(pose)  # Apply the trial mover for 1500 iterations
    mc.recover_low(pose)  # Recover the best pose after each trial move
    
    # Append the score to the list
    score = score3(pose)
    score_list.append(score)  # Store the score in the list
    if i % 1000 == 0:
        print(f"Iteration {i} - Score: {score}")  # Print the score after every 1000 iterations

# Recover best structure
mc.recover_low(pose)

print(score3(pose))

# Convert back to fullatom
to_fullatom = rosetta.protocols.simple_moves.SwitchResidueTypeSetMover("fa_standard")
to_fullatom.apply(pose)

# Output final structure
pose.dump_pdb("villin_abinitio.pdb")
print("Folding complete. Final structure saved as villin_abinitio.pdb")


## getting the native.pdb file
url = "https://files.rcsb.org/download/1VII.pdb"

response = requests.get(url)

with open("native.pdb", "w") as file:
    file.write(response.text)

print("Saved 1VII as native.pdb")

# plot the scores as a energy over iteration graph
import matplotlib.pyplot as plt

plt.plot(score_list)
plt.xlabel("Iteration")
plt.ylabel("Score")
plt.title("Energy vs Iteration")
plt.savefig("energy_vs_iteration.png")
# plt.show()
