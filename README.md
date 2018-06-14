# PyRAMA
Python3 implementation of the Ramachandran plot

*Package does work with Python2.7 but it is recommended to use Python3*

### For easy installation use pip:

    pip3 install pyrama
    

### Usage:

    pyrama my_pdb_file.pdb

Note: The script is able to read in multiple PDB files, however all of the torsion angles will be displayed on the same plot, with the same color!

### Example output:

![Example output](https://i.imgur.com/zOGxZ2r.png)

### Dependencies:

Running PyRAMA requires *matplotlib*, *numpy* and *biopython*

To install these on a standard Linux system:

    apt install python3-matplotlib
    apt install python3-biopython
    apt install python3-numpy

For the standard PSI and PHI preferences see:

Lovell *et al*. Structure validation by Calpha geometry: phi,psi and Cbeta deviation. 2003; DOI: 10.1002/prot.10286
