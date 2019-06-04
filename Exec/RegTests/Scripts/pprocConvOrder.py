#!/usr/bin/env python3

# Post-processing script for convergence analysis
# Must be used after multirun.py script
# Usage:
#   ./pprocConvOrder.py pproc_exec pproc_type

#   with pproc_exec the processing executable path:
#   - fcompare if pproc_type == 0. Analytical solution is known  
#   - diffsamedomain if pproc_type == 1. Analytical solution is not known and errors
#     are computed from the next finer grid  


import sys
import os
import fnmatch
import shutil
import numpy as np

import matplotlib
import matplotlib.pyplot as plt

def pproc(pproc_exe, pproc_type):

    # User data
    vars=["Y(CO2)", "y_velocity", "density", "Y(O2)", "Y(CH4)" ]
    resolution = [64,128,256,512]        

    # Get a local copy of post-processing executable
    run_dir = os.getcwd()
    if ( not os.path.isfile(os.path.basename(pproc_exe)) ):
        shutil.copy(pproc_exe, run_dir)
    test_name = run_dir.split("/")[-1]

    # Run the postprocessing
    if ( pproc_type == "0" ):     # running fcompare since analytical solution is known
        errors = np.empty([len(resolution),len(vars)+1])
        pltfile=[]
        for res in range(len(resolution)):
            case = resolution[res]
            errors[res,0] = case

            # Get the fcompare inputs: first and last solution of current case
            # TODO: the analytical solution might not be plt****_00000 ...
            for f in os.listdir(run_dir):
                if ( not fnmatch.fnmatch(f, '*old*')):
                    if (f.startswith("{}_plt_{}_".format(test_name,case))):
                        pltfile.append(f)
            pltfile.sort()
            outfile = "error_{}.analysis.out".format(case)
            os.system("./{} -n 2 {} {} > {}".format(os.path.basename(pproc_exe), pltfile[0], pltfile[-1], outfile))
            pltfile.clear()
        
            # Extract errors on each variable
            with open(outfile) as fp:
                for i, line in enumerate(fp):
                    if (i >= 5):
                        var = line.split()[0]
                        for v in range(len(vars)):
                            if ( var == vars[v] ):
                                errors[res,v+1] = line.split()[1]
            os.system("rm {}".format(outfile))
    elif ( pproc_type == "1" ):   # running diffsamedomain. No analytical sol ...
        errors = np.empty([len(resolution)-1,len(vars)+1])
        pltfile=[]
        pltfilenext=[]
        for res in range(len(resolution)-1):
            case = resolution[res]
            nextcase = resolution[res+1]
            errors[res,0] = case

            # Get the diffsamedomain inputs: last solutions of current 
            # and next finer cases. These run should have been runned to the same final time
            for f in os.listdir(run_dir):
                if ( not fnmatch.fnmatch(f, '*old*')):
                    if (f.startswith("{}_plt_{}_".format(test_name,case))):
                        pltfile.append(f)
                    if (f.startswith("{}_plt_{}_".format(test_name,nextcase))):
                        pltfilenext.append(f)
            pltfile.sort()
            pltfilenext.sort()
            outfile = "error_{}.analysis.out".format(case)
            os.system("./{} infile1={} reffile={} > {}".format(os.path.basename(pproc_exe), pltfile[-1], pltfilenext[-1], outfile))
            pltfile.clear()
            pltfilenext.clear()

            # Extract errors on each variable
            with open(outfile) as fp:
                for i, line in enumerate(fp):
                    if (i >= 5):
                        var = line.split(":")[0]
                        for v in range(len(vars)):
                            if ( var.split(" ")[0] == vars[v] ):
                                errors[res,v+1] = line.split(":")[1]
            os.system("rm {}".format(outfile))
    else:
        print("Wrong pproc_type: {}. should be either 0 or 1".format(pproc_type))
        return


    print(errors)
    # Plot data
    plotdata(errors, test_name, vars)
    writetex(errors, test_name, vars)

def plotdata(data, test_name, vars):
    # Evaluate 2nd order slope
    snd_order = data[:,1]*1.05
    for i in range(1,len(data[:,1])):
        snd_order[i] = snd_order[i-1]/np.exp(2.0*np.log(2.0))
    for i in range(0, len(vars)):    
        plt.plot(data[:,0], data[:,i+1], label="{}".format(vars[i]))
    plt.plot(data[:,0], snd_order[:],linestyle='--',color='k', label='2nd-order')
    plt.xlabel("Resolution")
    plt.ylabel("Error L2norm")
    plt.xscale("log")
    plt.yscale("log")
    plt.grid(which='both',color='k', linestyle=':', linewidth=1)
    plt.legend(bbox_to_anchor=(0.9, 0.9), loc=1, borderaxespad=0.)
    plt.savefig("Convergence_{}.png".format(test_name))

def writetex(data, test_name, vars):
    # Evaluate order
    conv_order = np.empty([len(data[:,0])-1,len(vars)])
    for v in range(len(vars)):
        for i in range(len(conv_order[:,0])):
            conv_order[i,v] = np.log(data[i,v+1]/data[i+1,v+1])/np.log(2.0)
    fout = open("ConvTable.tex", "w")            
    fout.write("\\begin{table}[ht!]\n")
    fout.write("\centering\n")
    fout.write("\\begin{tabular}{l|")
    for i in range(len(conv_order[:,0])):
        fout.write("c ")
    fout.write("}\n")    
    fout.write("\hline\n")
    fout.write("Variable ")
    for i in range(len(conv_order[:,0])):
        fout.write("&  {}/{} ".format(data[i+1,0],data[i,0]))
    fout.write("\\\\\n\hline\hline\n")
    for v in range(len(vars)):
        fout.write("{} ".format(vars[v].replace("_","\_")))
        for i in range(len(conv_order[:,0])):
            fout.write("&  {:.3f} ".format(conv_order[i,v]))
        fout.write("\\\\\n")
    fout.write("\end{tabular}\n")
    fout.write("\caption{PeleLM convergence order}\n")
    fout.write("\label{table:conv}\n")
    fout.write("\end{table}\n")
    fout.close()


if __name__ == "__main__":
    pproc_exe = str(sys.argv[1])
    pproc_type = str(sys.argv[2])
    pproc(pproc_exe, pproc_type)
