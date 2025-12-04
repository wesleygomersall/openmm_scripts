#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd

totalsteps = 25_000_000
time = 50
xaxis = 'Step'
yaxis = 'Distance(nm)'
xlabel = 'Time (ns)'
ylabel = 'Distance (nm)'
title = 'Displacement of peptide'

output_path = "allcomparison_distance.pdf" 

files = [["20250707T205103_ala_1/displacement.csv",
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
           "20250708T133217_c16m15_2/displacement.csv"]
         ]

colors = ["red", "blue", "green", "yellow", "black", "orange", "purple"]

labels = ["Alanines", "C16", "C16 mutant 12", "PMPNN negative", "C16 mutant 13", "C16 mutant 14", "C16 mutant 15"]

assert len(files) == len(colors) == len(labels)

for i, replicates in enumerate(files): 
    for j, file in enumerate(replicates):
        df = pd.read_csv(file)
        if j == 0: 
            plt.plot(time * df[xaxis] / totalsteps, df[yaxis], color=colors[i], label=labels[i])
        else: 
            plt.plot(time * df[xaxis] / totalsteps, df[yaxis], color=colors[i])

plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.title(title)
plt.legend()
plt.grid(False)
plt.savefig(output_path, format="pdf")
