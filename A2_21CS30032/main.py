import os

naccess_path = "for_student/"

# read the contents of this file and split the A and B chains in to separate pdb files
A_pdb_file = "A.pdb"
B_pdb_file = "B.pdb"
A_asa_file = "A.asa"
B_asa_file = "B.asa"

Complexes = {}

def select_sigma(A_pdb_line):
    l = A_pdb_line.split()
    sigma = 0
    if l[11] == 'N' or l[11] == 'O':
        sigma = -6
        if l[11] == 'N' and l[3] == "ARG" and l[2]!= 'N' :
            sigma = -50
        if l[11] == 'N' and l[3] == "LYS" and l[2]!= 'N' :
            sigma = -50
        if l[11] == 'O' and l[3] == "ASP" and l[2]!= 'O' :
            sigma = -24
        if l[11] == 'O' and l[3] == "GLU" and l[2]!= 'O' :
            sigma = -24

    elif l[11] == 'C':
        sigma = 16
    elif l[11] == 'S':
        sigma = 21

    return sigma


def compute_InterfaceAreaandSolvationEnergy(decoy_pdf,complex_no):
    # run naccess on the decoy pdb file
    os.system(naccess_path + "./naccess " + decoy_pdf)

    # compare the .asa files entries of decoy.asa with A.asa and B.asa
    # read the .asa files
    decoyfile = "complex"


    with open(decoyfile + ".asa", 'r') as f:
        decoy_lines = f.readlines()
    f.close()

    with open(A_asa_file, 'r') as f:
        A_lines = f.readlines()
    f.close()

    with open(B_asa_file, 'r') as f:
        B_lines = f.readlines()
    f.close()

    with open(A_pdb_file) as f :
        A_pdb_lines = f.readlines()
    f.close()

    with open(B_pdb_file) as f :
        B_pdb_lines = f.readlines()
    f.close()

    # compute the interface area and solvation_energy

    interface_area = 0
    solvation_energy = 0
    interface_check = []

    for i in range(len(decoy_lines)) :
        l = decoy_lines[i].split()
        if i < len(A_lines):
            diff = float(A_lines[i].split()[9]) - float(l[9])
            if diff > 0.1:
                interface_check.append(1)
                # interface area calculation
                interface_area += diff

                # solvation energy calculation
                sigma = select_sigma(A_pdb_lines[i])
                solvation_energy += diff * sigma
            else :
                interface_check.append(0)

        else:
            diff = float(B_lines[i-len(A_lines)].split()[9]) - float(l[9])
            if diff > 0.1:
                interface_check.append(1)
                # interface area calculation
                interface_area += diff

                # solvation energy calculation
                sigma = select_sigma(B_pdb_lines[i-len(A_lines)])
                solvation_energy += diff * sigma
            else :
                interface_check.append(0)

    # delete the files of the complex
    os.remove(decoyfile + ".asa")
    os.remove(decoyfile + ".rsa")
    os.remove(decoyfile + ".log")


    interface_area = interface_area / 2

    Complexes[complex_no] = {"Interface Area": interface_area}
    Complexes[complex_no]["Solvation Energy"] = solvation_energy
    Complexes[complex_no]["Interface Check"] = interface_check



target_pdb_file = "upload/target.pdb"



#create the files
with open(A_pdb_file, 'w') as f:
    f.write("")
f.close()
with open(B_pdb_file, 'w') as f:
    f.write("")
f.close()

with open(target_pdb_file, 'r') as f:
    lines = f.readlines()

    for line in lines :
        l = line.split()
        if l[4] == 'A':
            with open(A_pdb_file, 'a') as f:
                f.write(line)
            f.close()
        elif l[4] == 'B':
            with open(B_pdb_file, 'a') as f:
                f.write(line)
            f.close()

# run naccess on the A and B pdb files

os.system(naccess_path + "./naccess "+ A_pdb_file)
os.system(naccess_path + "./naccess " + B_pdb_file)


decoy_path = "upload/Decoys/"

for i in range(1,11) :
    complex_file = decoy_path + "complex." + str(i) + ".pdb"

    # compute interface area and solvation energy for each decoy
    compute_InterfaceAreaandSolvationEnergy(complex_file,i)


# determine the interface check array for the target pdb file
target_interface_check = []

# run naccess on the target pdb file
os.system(naccess_path + "./naccess " + target_pdb_file)

with open("target.asa", 'r') as f:
    target_lines = f.readlines()
f.close()

with open(A_asa_file, 'r') as f:
    A_lines = f.readlines()
f.close()

with open(B_asa_file, 'r') as f:
    B_lines = f.readlines()
f.close()

for i in range(len(target_lines)) :
    l = target_lines[i].split()
    if i < len(A_lines):
        diff = float(A_lines[i].split()[9]) - float(l[9])
        if diff > 0.1:
            target_interface_check.append(1)
        else :
            target_interface_check.append(0)

    else:
        diff = float(B_lines[i-len(A_lines)].split()[9]) - float(l[9])
        if diff > 0.1:
            target_interface_check.append(1)
        else :
            target_interface_check.append(0)





from Bio.PDB import PDBParser, Superimposer
import numpy as np

# computing LMRSD and IRMSD for all the decoys

for i in range(1,11) :
    parser = PDBParser(QUIET=True)
    decoy_path = f"upload/Decoys/complex.{i}.pdb"
    target_structure = parser.get_structure("target", "upload/target.pdb")
    decoy_structure = parser.get_structure("decoy", decoy_path)


    target_atoms_A = [atom for atom in target_structure[0]['A'].get_atoms()]
    decoy_atoms_A = [atom for atom in decoy_structure[0]['A'].get_atoms()]

    target_atoms_B = [atom for atom in target_structure[0]['B'].get_atoms()]
    decoy_atoms_B = [atom for atom in decoy_structure[0]['B'].get_atoms()]

    # Superimpose the receptor chain
    sup = Superimposer()
    sup.set_atoms(target_atoms_A, decoy_atoms_A)

    # Apply transformation to the receptor chain of decoy
    sup.apply(decoy_atoms_A)

    # apply the transformation to the ligand chain of the decoy
    sup.apply(decoy_atoms_B)

    # sup.set_atoms(target_atoms, decoy_atoms)

    # get the coordinates of the atoms of target and decoy
    target_coords_A = [atom.get_coord() for atom in target_atoms_A]
    target_coords_B = [atom.get_coord() for atom in target_atoms_B]
    decoy_coords_A = [atom.get_coord() for atom in decoy_atoms_A]
    decoy_coords_B = [atom.get_coord() for atom in decoy_atoms_B]

    target_coords = target_coords_A + target_coords_B
    decoy_coords = decoy_coords_A + decoy_coords_B

    # compute LRMSD and IRMSD and Fnat
    LRMSD = 0
    IRMSD = 0
    Fnat = 0

    for j in range(len(target_coords)):
        target_coord = target_coords[j]
        decoy_coord = decoy_coords[j]

        # add the sum of the square of the difference of the coordinates to the LRMSD
        LRMSD += np.sum((target_coord - decoy_coord) ** 2)

        if Complexes[i]["Interface Check"][j] == 1:
            IRMSD += np.sum((target_coord - decoy_coord) ** 2)
            if target_interface_check[j] == 1:
                Fnat += 1



    LRMSD = LRMSD / len(target_coords)
    IRMSD = IRMSD / len(target_coords)
    LRMSD = np.sqrt(LRMSD)
    IRMSD = np.sqrt(IRMSD)
    Fnat = Fnat / sum(target_interface_check)

    Complexes[i]["LRMSD"] = LRMSD
    Complexes[i]["IRMSD"] = IRMSD
    Complexes[i]["Fnat"] = Fnat

# for i in range(1,11) :
#     # print the lmrsd and irmsd for each decoy
#     print(f"Decoy {i}:")
#     print(f"    Interface Area: {Complexes[i]['Interface Area']}")
#     print(f"    Solvation Energy: {Complexes[i]['Solvation Energy']}")
#     print(f"    LRMSD: {Complexes[i]['LRMSD']}")
#     print(f"    IRMSD: {Complexes[i]['IRMSD']}")
#     print(f"    Fnat: {Complexes[i]['Fnat']}")


# tabulating the results
import tabulate
table = []
headers = ["Filename", "Interface Area", "Solvation Energy", "LRMSD", "IRMSD", "Fnat"]

for i in range(1,11) :
    table.append([f"complex.{i}.pdb", Complexes[i]["Interface Area"], Complexes[i]["Solvation Energy"], Complexes[i]["LRMSD"], Complexes[i]["IRMSD"], Complexes[i]["Fnat"]])

with open("Score.txt", 'w') as f:
    f.write(tabulate.tabulate(table, headers, tablefmt="grid"))

# cleaning up
os.remove(A_pdb_file)
os.remove(B_pdb_file)
os.remove("A.asa")
os.remove("B.asa")
os.remove("A.rsa")
os.remove("B.rsa")
os.remove("A.log")
os.remove("B.log")
os.remove("target.asa")
os.remove("target.rsa")
os.remove("target.log")
