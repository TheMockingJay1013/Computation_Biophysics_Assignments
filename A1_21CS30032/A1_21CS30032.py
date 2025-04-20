import requests
import re
import os

# a map of residue 3 letter code to 1 letter code
residue_map = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
    'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
    'TYR': 'Y', 'VAL': 'V'
}

AA_molecular_weight = {
    'ALA': 89.1,
    'ARG': 174.2,
    'ASN': 132.1,
    'ASP': 133.1,
    'CYS': 121.2,
    'GLU': 147.1,
    'GLN': 146.2,
    'GLY': 75.1,
    'HIS': 155.2,
    'ILE': 131.2,
    'LEU': 131.2,
    'LYS': 146.2,
    'MET': 149.2,
    'PHE': 165.2,
    'PRO': 115.1,
    'SER': 105.1,
    'THR': 119.1,
    'TRP': 204.2,
    'TYR': 181.2,
    'VAL': 117.1
}

chain_ASA = {}

def get_fasta_sequence(pdb_id):
    url_fasta = "https://www.rcsb.org/fasta/entry/" + pdb_id + "/display"
    fasta = requests.get(url_fasta)
    fasta_str = fasta.text
    sequences = {}
    order = []

    lines = fasta_str.strip().split("\n")
    current_chains = []

    for line in lines:
        if line.startswith(">"):
            match = re.search(r'Chains ([A-Z, ]+)\|', line)
            if match:
                chains = match.group(1).split(", ")
                current_chains = chains  # Store all chains linked to this sequence
                for chain in chains:
                    if chain not in sequences:
                        sequences[chain] = ""
                        order.append(chain)  # Maintain order of chains
        else:
            for chain in current_chains:
                sequences[chain] += line.strip()

    # Sort chains in proper order (A, B, C, D, etc.)
    order = sorted(order)

    # Concatenate sequences in chain order
    return "".join(sequences[chain] for chain in order if chain in sequences)

def get_pdb_structure(pdb_id):
    url_pdb = "https://files.rcsb.org/view/" + pdb_id + ".pdb"
    pdb = requests.get(url_pdb)

    structure = []
    for line in pdb.text.split("\n"):
        if line.startswith("ATOM"):
            l = line.split()
            d = {}
            d["atom"] = l[2]
            d['atom_num'] = int(l[1])
            d['residue'] = l[3]
            d['chain'] = l[4]
            d['residue_num'] = int(l[5])
            structure.append(d)

    return structure

def get_pdb_sequence(structure):
    pdb_seq = []
    l = [residue_map[structure[0]['residue']], structure[0]['chain']]
    pdb_seq.append(l)
    for i in range(1, len(structure)):
        if structure[i]['residue_num'] != structure[i-1]['residue_num']:
            l = []
            l.append(residue_map[structure[i]['residue']])
            l.append(structure[i]['chain'])
            pdb_seq.append(l)
    return pdb_seq

pdb_id = input("Enter the PDB ID: ")

# Open a text file to write the output
with open(f"A1_{pdb_id}.txt", "w") as f:
    # getting the sequence from the fasta response
    sequence = get_fasta_sequence(pdb_id)

    # getting the pdb structure
    pdb_structure = get_pdb_structure(pdb_id)

    # extract the sequence from the pdb structure
    pdb_seq = get_pdb_sequence(pdb_structure)

    # identify the number of chains in the structure
    chains = set()
    for d in pdb_structure:
        chains.add(d['chain'])

    f.write(f"The structure has {len(chains)} chains: {', '.join(chains)}\n")
    f.write("------------------------------------ xx ------------------------------------\n\n")

    
    # check for breaks in the chains
    f.write("Checking the pdb sequence for breaks...\n\n")
    i = 0
    j = 0
    f_flag = 1
    while i < len(sequence) and j < len(pdb_seq):
        if sequence[i] == pdb_seq[j][0]:
            i += 1
            j += 1
        else:
            f.write(f"{sequence[i]} missing in pdb sequence at position {i+1} in chain {pdb_seq[j][1]}\n")
            f_flag = 0
            i += 1

    if f_flag:
        f.write("All the residues in the sequence are present in the pdb structure\n")

    f.write("\n")
    f.write("------------------------------------ xx ------------------------------------\n\n")

    ## running NACCESS for computing ASA
    # load the pdb file 
    url_pdb = "https://files.rcsb.org/view/" + pdb_id + ".pdb"
    pdb = requests.get(url_pdb)
    pdb_file = open(f"{pdb_id}.pdb", "w")
    pdb_file.write(pdb.text)
    pdb_file.close()

    # add the path to the naccess folder in this variable
    naccess_folder = "for_student/"

    # run NACCESS 
    os.system(naccess_folder + f"./naccess {pdb_id}.pdb")

    for chain in chains:
        chain_ASA[chain] = {"ASA": 0, "Molecular Weight": 0, "Residues": 0}

    # read the rsa file
    rsa_file = open(f"{pdb_id}.rsa", "r")
    lines = rsa_file.readlines()
    for line in lines:
        if line.startswith("RES"):
            l = line.split()
            chain = l[2]
            chain_ASA[chain]["chain_id"] = chain
            chain_ASA[chain]["ASA"] += float(l[4])
            chain_ASA[chain]["Residues"] += 1
            chain_ASA[chain]["Molecular Weight"] += AA_molecular_weight[l[1]]

    # sort the chains based on the chain id
    chains = sorted(chains)
    
    import tabulate
    table = []
    # set the headers
    headers = ["Chain ID", "Residues", "Molecular Weight", "ASA"]
    for chain in chains:
        table.append([chain, chain_ASA[chain]["Residues"], chain_ASA[chain]["Molecular Weight"], chain_ASA[chain]["ASA"]])
        
    f.write(tabulate.tabulate(table, headers, tablefmt="grid"))
    f.write("\n")
    f.write("------------------------------------ xx ------------------------------------\n\n")
    
    # delete the pdb rsa asa and log files
    os.remove(f"{pdb_id}.rsa")
    os.remove(f"{pdb_id}.asa")
    os.remove(f"{pdb_id}.log")
    os.remove(f"{pdb_id}.pdb")
    