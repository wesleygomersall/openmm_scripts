#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd
import argparse

ydata_col = "tRNA_displacement(nm)"

myfiles = [["20250707T205103_ala_1/displacement.csv", 
            "20250707T205103_ala_2/displacement.csv", 
            "20250707T205103_ala_3/displacement.csv"], 
           ["20250707T205103_c16_1/displacement.csv", 
            "20250707T205103_c16_3/displacement.csv", 
            "20250707T205103_c16_2/displacement.csv"], 
           ["20250707T215351_c16m12_1/displacement.csv", 
            "20250707T215351_c16m12_2/displacement.csv", 
            "20250707T215351_c16m12_3/displacement.csv"], 
           ["20250707T215351_pmpnn_1/displacement.csv", 
            "20250707T215351_pmpnn_2/displacement.csv", 
            "20250707T215351_pmpnn_3/displacement.csv"], 
           ["20250708T133217_c16m13_2/displacement.csv", 
            "20250708T133217_c16m13_1/displacement.csv", 
            "20250708T133217_c16m13_3/displacement.csv"], 
           ["20250708T133217_c16m14_1/displacement.csv", 
            "20250708T133217_c16m14_2/displacement.csv", 
            "20250708T133217_c16m14_3/displacement.csv"], 
           ["20250708T133217_c16m15_1/displacement.csv", 
            "20250708T133217_c16m15_3/displacement.csv", 
            "20250708T133217_c16m15_2/displacement.csv"] ]
mycolors = ["red", "blue", "green", "yellow", "black", "orange", "purple"]
mylabels = ["Alanines", "C16", "C16 mutant 12", "PMPNN negative", "C16 mutant 13", "C16 mutant 14", "C16 mutant 15"]

def plot_dist(files, colors, labels, ydata, time = 50):   
    # time = 50 # in ns

    ydata = 'Distance(nm)' # TODO CHANGE!

    xlabel = 'Time (ns)'
    ylabel = 'Distance (nm)'
    title = 'Displacement of peptide'

    output_path = "allcomparison_distance.pdf" 

    assert len(files) == len(colors) == len(labels)

    for i, replicates in enumerate(files): 
        for j, file in enumerate(replicates):
            df = pd.read_csv(file)
            timestep = time / len(df.index)
            if j == 0: 
                plt.plot(timestep * df.index+1, df[ydata], color=colors[i], label=labels[i])
            else: 
                plt.plot(timestep * df.index+1, df[ydata], color=colors[i])

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.grid(False)
    plt.savefig(output_path, format="pdf")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--datacol", "-d", type=str, default = "colname", help="Column name with position data")
    parser.add_argument("--time", "-t", type=float, default = 50.0, help="Time in nanoseconds")
    args = parser.parse_args() 

    args.datacol = ydata_col

    plot_dist(myfiles, mycolors, mylabels, args.datacol, args.time)
