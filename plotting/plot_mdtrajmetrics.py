#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd
import argparse

mycolors = ["red", "blue"]
mylabels = ["M8", "M8_14"]
mystyles = ["solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid", "solid"]

displacement_col = "Ligand Displacement" # column titles in csvs
# rmsd_col = "Peptide_RMSD(A)"

disp_1 = [
"20260122T095757_comdn230_M8_rep2/output_displacement.csv",
"20260122T095757_comdn230_M8_rep3/output_displacement.csv",
"20260122T095805_comdn230_M8_rep1/output_displacement.csv",
        ] 
disp_2 = [
"20260122T095757_comdn230_M8-14_rep1/output_displacement.csv",
"20260122T095757_comdn230_M8-14_rep2/output_displacement.csv",
"20260122T095757_comdn230_M8-14_rep3/output_displacement.csv",
        ] 

rmsd_1 = [
"20260122T095757_comdn230_M8_rep2/output_RMSDs.csv",
"20260122T095757_comdn230_M8_rep3/output_RMSDs.csv",
"20260122T095805_comdn230_M8_rep1/output_RMSDs.csv",
        ] 
rmsd_2 = [
"20260122T095757_comdn230_M8-14_rep1/output_RMSDs.csv",
"20260122T095757_comdn230_M8-14_rep2/output_RMSDs.csv",
"20260122T095757_comdn230_M8-14_rep3/output_RMSDs.csv",
        ] 

def plot_mdtraj_outputs(files, colors, labels, ydata, time, 
                        ylabel = '',
                        title = '',
                        output_path = 'mdtraj_output_UNKNOWN.pdf',
                        xlabel = 'Time (ns)',
                        linestyles = ["solid"]):   

    plt.clf() 

    if not output_path.endswith('.pdf'):
        output_path = output_path + '.pdf'

    assert len(files) == len(colors) == len(labels)

    for i, replicates in enumerate(files): 
        for j, file in enumerate(replicates):
            df = pd.read_csv(file)
            timestep = time / len(df.index)
            if j == 0: 
                plt.plot(timestep * df.index+1, df[ydata], color=colors[i], linestyle=linestyles[0], label=labels[i])
            else: 
                try: # get a valid linestyle 
                    linestyles[j]
                except IndexError: 
                    print("linestyles don't match groupings (there will be duplicated colors and styles)")
                    style = linestyles[-1]
                else:
                    style = linestyles[j]

                plt.plot(timestep * df.index+1, df[ydata], color=colors[i], linestyle=style)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.grid(False)
    plt.savefig(output_path, format="pdf")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--datacol", "-d", type=str, default = "Chain0_displacement(nm)", help="Column name for column of displacement data per trajectory frame. Default 'Chain0_displacement(nm)'.")
    parser.add_argument("--time", "-t", type=float, default = 50.0, help="Time of trajectory in nanoseconds, default 50 ns.")
    args = parser.parse_args() 

    all_displacements = []
    all_displacements.append(disp_1)
    all_displacements.append(disp_2) 

    plot_mdtraj_outputs(all_displacements,
                        mycolors, mylabels, displacement_col, args.time, 
                        ylabel = 'Displacement (nm)',
                        title = 'Peptide displacement',
                        output_path = 'mdtraj_output_displacement.pdf',
                        linestyles = mystyles)

    '''
    all_rmsds = []
    all_rmsds.append(rmsd_1)
    all_rmsds.append(rmsd_2) 
    args.datacol = rmsd_col
    plot_mdtraj_outputs(all_rmsds, 
                        mycolors, mylabels, rmsd_col, args.time, 
                        ylabel = 'RMSD (A)',
                        title = 'Peptide RMSD',
                        output_path = 'mdtraj_output_rmsd.pdf',
                        linestyles = mystyles)

    '''
