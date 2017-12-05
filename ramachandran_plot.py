import math
import sys
import os

import matplotlib.pyplot as plt
import numpy as np
from Bio import PDB
from matplotlib import colors

"""
Written by Gábor Erdős, 2017
Contact info: gerdos[at]caesar.elte.hu

The preferences were calculated from the following artice:
Lovell et al. Structure validation by Calpha geometry: phi,psi and Cbeta deviation. 2003
DOI: 10.1002/prot.10286
"""

rama_preferences = {
    "General": {
        "file": "pref_general.data",
        "cmap": colors.ListedColormap(['#FFFFFF', '#B3E8FF', '#7FD9FF']),
        "bounds": [0, 0.0005, 0.02, 1],
    },
    "GLY": {
        "file": "pref_glycine.data",
        "cmap": colors.ListedColormap(['#FFFFFF', '#FFE8C5', '#FFCC7F']),
        "bounds": [0, 0.002, 0.02, 1],
    },
    "PRO": {
        "file": "pref_proline.data",
        "cmap": colors.ListedColormap(['#FFFFFF', '#D0FFC5', '#7FFF8C']),
        "bounds": [0, 0.002, 0.02, 1],
    },
    "PRE-PRO": {
        "file": "pref_preproline.data",
        "cmap": colors.ListedColormap(['#FFFFFF', '#B3E8FF', '#7FD9FF']),
        "bounds": [0, 0.002, 0.02, 1],
    }
}
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
rama_pref_values = {}
for key, val in rama_preferences.items():
    rama_pref_values[key] = np.full((360, 360), 0, dtype=np.float64)
    with open(os.path.join(__location__, val["file"])) as fn:
        for line in fn:
            if not line.startswith("#"):
                # Preference file has values for every second position only
                rama_pref_values[key][int(float(line.split()[1])) + 180][int(float(line.split()[0])) + 180] = float(
                    line.split()[2])
                rama_pref_values[key][int(float(line.split()[1])) + 179][int(float(line.split()[0])) + 179] = float(
                    line.split()[2])
                rama_pref_values[key][int(float(line.split()[1])) + 179][int(float(line.split()[0])) + 180] = float(
                    line.split()[2])
                rama_pref_values[key][int(float(line.split()[1])) + 180][int(float(line.split()[0])) + 179] = float(
                    line.split()[2])

normals = {}
outliers = {}
for key, val in rama_preferences.items():
    normals[key] = {"x": [], "y": []}
    outliers[key] = {"x": [], "y": []}
for inp in sys.argv[1:]:
    if not os.path.isfile(inp):
        print("{} not found!".format(inp))
        continue
    structure = PDB.PDBParser().get_structure('input_structure', inp)
    for model in structure:
        for chain in model:
            polypeptides = PDB.PPBuilder().build_peptides(chain)
            for poly_index, poly in enumerate(polypeptides):
                phi_psi = poly.get_phi_psi_list()
                for res_index, residue in enumerate(poly):
                    res_name = "{}".format(residue.resname)
                    res_num = residue.id[1]
                    phi, psi = phi_psi[res_index]
                    if phi and psi:
                        aa_type = ""
                        if str(poly[res_index + 1].resname) == "PRO":
                            aa_type = "PRE-PRO"
                        elif res_name == "PRO":
                            aa_type = "PRO"
                        elif res_name == "GLY":
                            aa_type = "GLY"
                        else:
                            aa_type = "General"
                        if rama_pref_values[aa_type][int(math.degrees(psi)) + 180][int(math.degrees(phi)) + 180] < \
                                rama_preferences[aa_type]["bounds"][1]:
                            print("{} {} {} {}{} is an outlier".format(inp, model, chain, res_name, res_num))
                            outliers[aa_type]["x"].append(math.degrees(phi))
                            outliers[aa_type]["y"].append(math.degrees(psi))
                        else:
                            normals[aa_type]["x"].append(math.degrees(phi))
                            normals[aa_type]["y"].append(math.degrees(psi))

for idx, (key, val) in enumerate(sorted(rama_preferences.items(), key=lambda x: x[0].lower())):
    plt.subplot(2, 2, idx + 1)
    plt.title(key)
    plt.imshow(rama_pref_values[key], cmap=rama_preferences[key]["cmap"],
               norm=colors.BoundaryNorm(rama_preferences[key]["bounds"], rama_preferences[key]["cmap"].N),
               extent=(-180, 180, 180, -180))
    plt.scatter(normals[key]["x"], normals[key]["y"])
    plt.scatter(outliers[key]["x"], outliers[key]["y"], color="red")
    plt.xlim([-180, 180])
    plt.ylim([-180, 180])
    plt.plot([-180, 180], [0, 0], color="black")
    plt.plot([0, 0], [-180, 180], color="black")
    plt.locator_params(axis='x', nbins=7)
    plt.xlabel(r'$\phi$')
    plt.ylabel(r'$\psi$')
    plt.grid()

plt.tight_layout()
# plt.savefig("asd.png", dpi=300)
plt.show()
