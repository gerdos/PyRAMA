from __future__ import division, print_function

import math
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from Bio import PDB
from matplotlib import colors

from .config import RAMA_PREFERENCES

RAMA_PREF_VALUES = None


def _cache_RAMA_PREF_VALUES():
    f_path = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
    RAMA_PREF_VALUES = {}
    for key, val in RAMA_PREFERENCES.items():
        RAMA_PREF_VALUES[key] = np.full((360, 360), 0, dtype=np.float64)
        with open(os.path.join(f_path, val["file"])) as fn:
            for line in fn:
                if line.startswith("#"):
                    continue
                else:
                    x = int(float(line.split()[1]))
                    y = int(float(line.split()[0]))
                    RAMA_PREF_VALUES[key][x + 180][y + 180] \
                        = RAMA_PREF_VALUES[key][x + 179][y + 179] \
                        = RAMA_PREF_VALUES[key][x + 179][y + 180] \
                        = RAMA_PREF_VALUES[key][x + 180][y + 179] \
                        = float(line.split()[2]) 
    return RAMA_PREF_VALUES


def calc_ramachandran(file_name_list):
    """
    Main calculation and plotting definition
    :param file_name_list: List of PDB files to plot
    :return: results dict by chain, with 'normals'/'outliers' keys for each chain
    """
    global RAMA_PREF_VALUES

    if RAMA_PREF_VALUES is None:
        RAMA_PREF_VALUES = _cache_RAMA_PREF_VALUES()

    results = {}
    # Calculate the torsion angle of the inputs
    for inp in file_name_list:
        if not os.path.isfile(inp):
            continue
        structure = PDB.PDBParser().get_structure('input_structure', inp)
        for model in structure:
            for chain in model:
                chain_id = chain.id
                # 初始化结构，兼容plot_ramachandran
                if chain_id not in results:
                    results[chain_id] = {"normals": {}, "outliers": {}}
                    for key in RAMA_PREFERENCES.keys():
                        results[chain_id]["normals"][key] = {"x": [], "y": [], "residues": []}
                        results[chain_id]["outliers"][key] = {"x": [], "y": [], "residues": []}
                polypeptides = PDB.PPBuilder().build_peptides(chain)
                for poly_index, poly in enumerate(polypeptides):
                    phi_psi = poly.get_phi_psi_list()
                    for res_index, residue in enumerate(poly):
                        res_name = "{}".format(residue.resname)
                        res_num = residue.id[1]
                        phi, psi = phi_psi[res_index]
                        if phi and psi:
                            if res_index + 1 < len(poly) and str(poly[res_index + 1].resname) == "PRO":
                                aa_type = "PRE-PRO"
                            elif res_name == "PRO":
                                aa_type = "PRO"
                            elif res_name == "GLY":
                                aa_type = "GLY"
                            else:
                                aa_type = "General"
                            phi_deg = math.degrees(phi)
                            psi_deg = math.degrees(psi)
                            residue_info = {"name": res_name, "number": res_num}
                            if RAMA_PREF_VALUES[aa_type][int(psi_deg) + 180][int(phi_deg) + 180] < \
                                    RAMA_PREFERENCES[aa_type]["bounds"][1]:
                                results[chain_id]["outliers"][aa_type]["x"].append(phi_deg)
                                results[chain_id]["outliers"][aa_type]["y"].append(psi_deg)
                                results[chain_id]["outliers"][aa_type]["residues"].append(residue_info)
                            else:
                                results[chain_id]["normals"][aa_type]["x"].append(phi_deg)
                                results[chain_id]["normals"][aa_type]["y"].append(psi_deg)
                                results[chain_id]["normals"][aa_type]["residues"].append(residue_info)
    return results


def plot_ramachandran(results, output_dir=None, show_plots=True):
    """
    Plot Ramachandran plots for each chain separately
    :param results: Dictionary with chain-specific data
    :param output_dir: Directory to save plot files (if None, won't save)
    :param show_plots: Whether to display plots (True) or just save them (False)
    :return: List of paths to saved plot files
    """
    global RAMA_PREF_VALUES
    if RAMA_PREF_VALUES is None:
        RAMA_PREF_VALUES = _cache_RAMA_PREF_VALUES()
    
    saved_files = []
    
    # Create a separate figure for each chain
    for chain_id, chain_data in results.items():
        plt.figure(figsize=(12, 10))
        plt.suptitle(f"Ramachandran Plot - Chain {chain_id}")
        
        normals = chain_data["normals"]
        outliers = chain_data["outliers"]
        
        for idx, (key, val) in enumerate(sorted(RAMA_PREFERENCES.items(), key=lambda x: x[0].lower())):
            plt.subplot(2, 2, idx + 1)
            plt.title(key)
            plt.imshow(RAMA_PREF_VALUES[key], cmap=RAMA_PREFERENCES[key]["cmap"],
                    norm=colors.BoundaryNorm(RAMA_PREFERENCES[key]["bounds"], RAMA_PREFERENCES[key]["cmap"].N),
                    extent=(-180, 180, 180, -180))
            
            # Plot normal points
            if normals[key]["x"]:
                plt.scatter(normals[key]["x"], normals[key]["y"], label="Normal")
            
            # Plot outliers
            if outliers[key]["x"]:
                plt.scatter(outliers[key]["x"], outliers[key]["y"], color="red", label="Outlier")
                # Add annotations for outliers
                for x, y, res in zip(outliers[key]["x"], outliers[key]["y"], outliers[key]["residues"]):
                    plt.annotate(f"{res['name']}{res['number']}", (x, y), textcoords="offset points", xytext=(0, 5), ha='center', fontsize=7, color='red')
            
            plt.xlim([-180, 180])
            plt.ylim([-180, 180])
            plt.plot([-180, 180], [0, 0], color="black")
            plt.plot([0, 0], [-180, 180], color="black")
            plt.locator_params(axis='x', nbins=7)
            plt.xlabel(r'$\phi$')
            plt.ylabel(r'$\psi$')
            plt.grid()
            
            if normals[key]["x"] or outliers[key]["x"]:
                plt.legend()

        plt.tight_layout()
        
        # Save plot if output directory is provided
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
            plot_file = os.path.join(output_dir, f"ramachandran_chain_{chain_id}.png")
            plt.savefig(plot_file, dpi=300)
            saved_files.append(plot_file)
        
        if show_plots:
            plt.show()
        else:
            plt.close()
    
    return saved_files
