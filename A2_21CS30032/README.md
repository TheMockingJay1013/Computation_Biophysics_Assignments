
The main.py is the python script to be run.
Before running it, ensure the following :
  1. The 'naccess_path' variable  (can be found and the beginning of the main.py file) is set to the correct folder containing the naccess
  2. The 'upload' folder given in the assignment is also included in the same folder that main.py is in.
  3. Biopython, numpy and tabulate libraries are installed in the system/environment.

The code will create several files as it executes, however the code to remove these files are also included in main.py .
If you wish to see these files, simply comment out the os.remove() of that file at the end of main.py


### Implementation details
For getting the delta ASA for the target of the decoys, we needed a reference ASA of the individual chains.
To obtain this, the code splits the chains A and B in the target.pdb file and creates 2 new pdb files A.pdb and B.pdb
It runs naccess on these files to get the ASA of the chains.
The code then runs naccess on the decoys and computes the change with the above outputs to compute the interface area.
An atom is considered an interface atom if its delta ASA is greater than 0.1

For finding the RMSD values,
First the receptors chains of the target and decoy are aligned using the Biopython Superimposer() object.
Then the receptor and ligand of the decoy is aligned to the target complex by applying the rotation and translation obtained from Superimposer()
LRMSD and IRMSD are then computed on the aligned complexes.

For finding the Fnat scores
For every decoy complex and the target complex, an array is created called the 'interface_check' (made during the interface area calculation) which is a binary array that tells if an atom is interface atom or not
This array is used to compute the fnat score
