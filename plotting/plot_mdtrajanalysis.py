#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd

totalsteps = 25_000_000
time = 50 # in nanoseconds
frames = 100 # 100 frames of pdb


c16ala_files = [["20250707T205103_ala_1/mdtrajanalysis.csv",
                 "20250707T205103_ala_2/mdtrajanalysis.csv",
                 "20250707T205103_ala_3/mdtrajanalysis.csv"],
                ["20250707T205103_c16_1/mdtrajanalysis.csv",
                 "20250707T205103_c16_2/mdtrajanalysis.csv",
                 "20250707T205103_c16_3/mdtrajanalysis.csv"],
                ["20250707T215351_pmpnn_1/mdtrajanalysis.csv",
                 "20250707T215351_pmpnn_2/mdtrajanalysis.csv",
                 "20250707T215351_pmpnn_3/mdtrajanalysis.csv"]]
c16ala_colors = ["red", "blue", "yellow"]
c16ala_labels = ["Alanines", "C16", "PMPNN"]
assert len(c16ala_files) == len(c16ala_colors) == len(c16ala_labels)


myfiles = [["20250707T215351_c16m12_1/mdtrajanalysis.csv",
            "20250707T215351_c16m12_2/mdtrajanalysis.csv",
            "20250707T215351_c16m12_3/mdtrajanalysis.csv"],
           ["20250708T133217_c16m13_1/mdtrajanalysis.csv",
            "20250708T133217_c16m13_2/mdtrajanalysis.csv",
            "20250708T133217_c16m13_3/mdtrajanalysis.csv"],
           ["20250708T133217_c16m14_1/mdtrajanalysis.csv",
            "20250708T133217_c16m14_2/mdtrajanalysis.csv",
            "20250708T133217_c16m14_3/mdtrajanalysis.csv"],
           ["20250708T133217_c16m15_1/mdtrajanalysis.csv",
            "20250708T133217_c16m15_2/mdtrajanalysis.csv",
            "20250708T133217_c16m15_3/mdtrajanalysis.csv"]]
mycolors = ["green", "black", "orange", "purple"]
mylabels = ["C16 mutant 12", "C16 mutant 13", "C16 mutant 14", "C16 mutant 15"]
assert len(myfiles) == len(mycolors) == len(mylabels)


def plot_simulation(files: list, colors: list, labels: list, 
                    graph_type: str, optional_label: str = ""): 
    plt.clf() # clear figure from any previous plotting
    xaxis = 'Frame'
    xlabel = 'Time (ns)'
    # valid_graph_types = ["displacement", "rmsd", "a-carbons"]
    match graph_type: 
        case "displacement":
            yaxis = 'Displacement(nm)'
            ylabel = 'Distance (nm)'
            title = 'Displacement of peptide'
        case "rmsd":
            yaxis = 'RMSD(nm)'
            ylabel = 'RMSD (nm)'
            title = 'Peptide RMSD'
        case "a-carbons": 
            yaxis = 'aCarbon-aCarbon Distance(nm)'
            ylabel = 'Distance (nm)'
            title = 'aCarbon-aCarbon Distance'
        case _: 
            raise ValueError("Use a valid graph_type")

    for i, replicates in enumerate(files):
        for j, file in enumerate(replicates):
            df = pd.read_csv(file)
            if j == 0: 
                plt.plot(time * df[xaxis] / frames, df[yaxis], color=colors[i], label=labels[i])
            else: 
                plt.plot(time * df[xaxis] / frames, df[yaxis], color=colors[i])

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.grid(False)
    if optional_label != "": optional_label = optional_label + "_"
    output_path = f"{optional_label}comparison_{graph_type}.pdf" 
    plt.savefig(output_path, format="pdf")

def plot_all(files: list, colors: list, labels: list,  optional_label: str = ""): 
    plot_simulation(files, colors, labels, "displacement", optional_label)
    plot_simulation(files, colors, labels, "rmsd", optional_label)
    plot_simulation(files, colors, labels, "a-carbons", optional_label)

if __name__ == "__main__":
    plot_all(myfiles, mycolors, mylabels, "mutants")
    plot_all(c16ala_files, c16ala_colors, c16ala_labels, "c16ala")
