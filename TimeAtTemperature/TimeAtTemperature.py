
import sys
import numpy as np     # installed with matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from scipy.integrate import cumtrapz
from mpl_toolkits.mplot3d import Axes3D

def main():
    # Read in datafiles
    ds1 = pd.read_csv('../data/sample - 1deg per min, 7.20mg.txt', sep='\s+', skiprows=lambda x: x == 1)
    ds2 = pd.read_csv('../data/sample - 2deg per min, 6.82mg.txt', sep='\s+', skiprows=lambda x: x == 1)
    ds3 = pd.read_csv('../data/sample - 5deg per min, 11.50mg.txt', sep='\s+', skiprows=lambda x: x == 1)
    ds4 = pd.read_csv('../data/sample - 7deg per min, 7.20mg.txt', sep='\s+', skiprows=lambda x: x == 1)
    ds5 = pd.read_csv('../data/sample - 30deg per min, 7.00mg.txt', sep='\s+', skiprows=lambda x: x == 1)

    # Correct "value" field (Heat flux) to (Heat flux per unit mass)
    ds1['Value'] /= 7.2
    ds2['Value'] /= 6.82
    ds3['Value'] /= 11.5
    ds4['Value'] /= 7.2
    ds5['Value'] /= 7.0

    plt.ion()

    # Plot heat flux vs time
    plt.figure(1)
    plt.plot(ds1['t'], ds1['Value'])
    plt.plot(ds2['t'], ds2['Value'])
    plt.plot(ds3['t'], ds3['Value'])
    plt.plot(ds4['t'], ds4['Value'])
    plt.plot(ds5['t'], ds5['Value'])
    plt.show()

    # Calculate total heat of reaction by integrating heat flux over time
    heat1 = cumtrapz(ds1['Value'],ds1['t'],initial = 0)
    heat2 = cumtrapz(ds2['Value'],ds2['t'],initial = 0)
    heat3 = cumtrapz(ds3['Value'],ds3['t'],initial = 0)
    heat4 = cumtrapz(ds4['Value'],ds4['t'],initial = 0)
    heat5 = cumtrapz(ds5['Value'],ds5['t'],initial = 0)

    # Plot total heat of reaction
    plt.figure(2)
    plt.plot(ds1['t'], heat1)
    plt.plot(ds2['t'], heat2)
    plt.plot(ds3['t'], heat3)
    plt.plot(ds4['t'], heat4)
    plt.plot(ds5['t'], heat5)
    plt.show()

    # Calculate degree of cure
    degree_of_cure1 = heat1 / heat1[-1]
    degree_of_cure2 = heat2 / heat2[-1]
    degree_of_cure3 = heat3 / heat3[-1]
    degree_of_cure4 = heat4 / heat4[-1]
    degree_of_cure5 = heat5 / heat5[-1]

    # Plot degree of cure
    plt.figure(3)
    plt.plot(ds1['t'], degree_of_cure1)
    plt.plot(ds2['t'], degree_of_cure2)
    plt.plot(ds3['t'], degree_of_cure3)
    plt.plot(ds4['t'], degree_of_cure4)
    plt.plot(ds5['t'], degree_of_cure5)
    plt.show()

    # Calculate rate of cure, d(alpha)/dt
    rate_of_cure1 = np.zeros(degree_of_cure1.shape)
    rate_of_cure1[:-1] = np.diff(degree_of_cure1) / np.diff(ds1['t'])
    rate_of_cure1[-1] = rate_of_cure1[-2]
    
    rate_of_cure2 = np.zeros(degree_of_cure2.shape)
    rate_of_cure2[:-1] = np.diff(degree_of_cure2) / np.diff(ds2['t'])
    rate_of_cure2[-1] = rate_of_cure2[-2]

    rate_of_cure3 = np.zeros(degree_of_cure3.shape)
    rate_of_cure3[:-1] = np.diff(degree_of_cure3) / np.diff(ds3['t'])
    rate_of_cure3[-1] = rate_of_cure3[-2]

    rate_of_cure4 = np.zeros(degree_of_cure4.shape)
    rate_of_cure4[:-1] = np.diff(degree_of_cure4) / np.diff(ds4['t'])
    rate_of_cure4[-1] = rate_of_cure4[-2]

    rate_of_cure5 = np.zeros(degree_of_cure5.shape)
    rate_of_cure5[:-1] = np.diff(degree_of_cure5) / np.diff(ds5['t'])
    rate_of_cure5[-1] = rate_of_cure5[-2]

    plt.ioff()

    # Plot rate of cure vs degree of cure AND temperature
    fig = plt.figure(4)
    ax = fig.add_subplot(111,projection='3d')
    ax.scatter(ds1['Ts'], degree_of_cure1, rate_of_cure1)
    ax.scatter(ds2['Ts'], degree_of_cure2, rate_of_cure2)
    ax.scatter(ds3['Ts'], degree_of_cure3, rate_of_cure3)
    ax.scatter(ds4['Ts'], degree_of_cure4, rate_of_cure4)
    ax.scatter(ds5['Ts'], degree_of_cure5, rate_of_cure5)
    plt.show()

main()