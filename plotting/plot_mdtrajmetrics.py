#!/usr/bin/env python3

import matplotlib.pyplot as plt
import pandas as pd
import argparse

mycolors = ["red", "blue", "green", "yellow", "orange", "purple", "black"]
mylabels = ["0.1", "0.25", "0.5", "1", "2", "5", "10"]
mystyles = ["dotted", "solid"]
displacement_col = "tRNA_displacement(nm)"
displacement_files = [["/Users/wesg/PHlab/data/20251201_pus4forcetests_SMD/20251201T141527_mut11_forcetest_0.1/output_displacement.csv", 
                       "/Users/wesg/PHlab/data/20251201_pus4forcetests_SMD/20251201T141527_wt_forcetest_0.1/output_displacement.csv",], 
                      ["/Users/wesg/PHlab/data/20251201_pus4forcetests_SMD/20251201T141527_mut11_forcetest_0.25/output_displacement.csv", 
                       "/Users/wesg/PHlab/data/20251201_pus4forcetests_SMD/20251201T141527_wt_forcetest_0.25/output_displacement.csv",], 
                      ["/Users/wesg/PHlab/data/20251201_pus4forcetests_SMD/20251201T141527_mut11_forcetest_0.5/output_displacement.csv", 
                       "/Users/wesg/PHlab/data/20251201_pus4forcetests_SMD/20251201T141527_wt_forcetest_0.5/output_displacement.csv",], 
                      ["/Users/wesg/PHlab/data/20251201_pus4forcetests_SMD/20251201T141527_mut11_forcetest_1/output_displacement.csv", 
                       "/Users/wesg/PHlab/data/20251201_pus4forcetests_SMD/20251201T141527_wt_forcetest_1/output_displacement.csv",],
                      ["/Users/wesg/PHlab/data/20251201_pus4forcetests_SMD/20251201T141527_mut11_forcetest_2/output_displacement.csv", 
                       "/Users/wesg/PHlab/data/20251201_pus4forcetests_SMD/20251201T141527_wt_forcetest_2/output_displacement.csv",],
                      ["/Users/wesg/PHlab/data/20251201_pus4forcetests_SMD/20251201T141527_mut11_forcetest_5/output_displacement.csv", 
                       "/Users/wesg/PHlab/data/20251201_pus4forcetests_SMD/20251201T141527_wt_forcetest_5/output_displacement.csv",],
                      ["/Users/wesg/PHlab/data/20251201_pus4forcetests_SMD/20251201T141528_mut11_forcetest_10/output_displacement.csv", 
                       "/Users/wesg/PHlab/data/20251201_pus4forcetests_SMD/20251201T141528_wt_forcetest_10/output_displacement.csv"],
                      ]

rmsd_col = "tRNA_RMSD(A)"
rmsd_files = [["/Users/wesg/PHlab/data/20251201_pus4forcetests_SMD/20251201T141527_mut11_forcetest_0.1/output_RMSDs.csv", 
               "/Users/wesg/PHlab/data/20251201_pus4forcetests_SMD/20251201T141527_wt_forcetest_0.1/output_RMSDs.csv"], 
              ["/Users/wesg/PHlab/data/20251201_pus4forcetests_SMD/20251201T141527_mut11_forcetest_0.25/output_RMSDs.csv", 
               "/Users/wesg/PHlab/data/20251201_pus4forcetests_SMD/20251201T141527_wt_forcetest_0.25/output_RMSDs.csv"], 
              ["/Users/wesg/PHlab/data/20251201_pus4forcetests_SMD/20251201T141527_mut11_forcetest_0.5/output_RMSDs.csv", 
               "/Users/wesg/PHlab/data/20251201_pus4forcetests_SMD/20251201T141527_wt_forcetest_0.5/output_RMSDs.csv"], 
              ["/Users/wesg/PHlab/data/20251201_pus4forcetests_SMD/20251201T141527_mut11_forcetest_1/output_RMSDs.csv", 
               "/Users/wesg/PHlab/data/20251201_pus4forcetests_SMD/20251201T141527_wt_forcetest_1/output_RMSDs.csv"], 
              ["/Users/wesg/PHlab/data/20251201_pus4forcetests_SMD/20251201T141527_mut11_forcetest_2/output_RMSDs.csv", 
               "/Users/wesg/PHlab/data/20251201_pus4forcetests_SMD/20251201T141527_wt_forcetest_2/output_RMSDs.csv"], 
              ["/Users/wesg/PHlab/data/20251201_pus4forcetests_SMD/20251201T141527_mut11_forcetest_5/output_RMSDs.csv", 
               "/Users/wesg/PHlab/data/20251201_pus4forcetests_SMD/20251201T141527_wt_forcetest_5/output_RMSDs.csv"], 
              ["/Users/wesg/PHlab/data/20251201_pus4forcetests_SMD/20251201T141528_mut11_forcetest_10/output_RMSDs.csv" , 
               "/Users/wesg/PHlab/data/20251201_pus4forcetests_SMD/20251201T141528_wt_forcetest_10/output_RMSDs.csv"], 
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

    args.datacol = displacement_col

    plot_mdtraj_outputs(displacement_files, 
                        mycolors, mylabels, args.datacol, args.time, 
                        ylabel = 'Displacement (nm)',
                        title = 'Pus4 force tests tRNA displacement',
                        output_path = 'mdtraj_output_displacement.pdf',
                        linestyles = mystyles)

    args.datacol = rmsd_col
    plot_mdtraj_outputs(rmsd_files, 
                        mycolors, mylabels, args.datacol, args.time, 
                        ylabel = 'RMSD (A)',
                        title = 'Pus4 force tests tRNA RMSD',
                        output_path = 'mdtraj_output_rmsd.pdf',
                        linestyles = mystyles)
