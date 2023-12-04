import sys
import argparse 
import matplotlib.pyplot as plt
import numpy as np

#plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['text.usetex'] = True
plt.rcParams['axes.formatter.use_mathtext'] = True
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['figure.autolayout'] = True
plt.rcParams['grid.linewidth'] = 0.5
figsize = plt.rcParams['figure.figsize']
figsize[0] *= 1.4
plt.rcParams['figure.figsize'] = figsize
fs = 'x-large'
plt.rcParams['axes.labelsize'] = fs
plt.rcParams['legend.fontsize'] = fs
plt.rcParams['xtick.labelsize'] = fs
plt.rcParams['ytick.labelsize'] = fs

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('file_csv', type=str, help='a .csv file to plot')
    args = parser.parse_args()

    try:
        data = np.loadtxt(args.file_csv, delimiter=',')
    except:
        print('Error: Cannot open csv file')
        quit()

    if (args.file_csv.strip()[-6:-4] == 'f1'):
        plt.plot(data[:,0], data[:,1], color='b', label=r'$\bf{mm\_kij}$')
        plt.plot(data[:,0], data[:,2], color='r', label=r'$\bf{mm\_jki}$')
    else:
        plt.plot(data[:,0], data[:,1], color='b', label=r'$\bf{mm\_blocks\_a}$')
        plt.plot(data[:,0], data[:,2], color='r', label=r'$\bf{mm\_blocks\_b}$')
        
    plt.xlabel(r'$N$')
    plt.ylabel(r'$\mathrm{MFLOPS}$')
    plt.grid()
    plt.legend(loc='best')
    plt.savefig(args.file_csv.strip()[:-3] + 'png')