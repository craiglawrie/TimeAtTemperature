"""Calculate "Time At Temperature" required for a thermosetting resin to reach a desired degree of cure.
   Takes as input data from Differential Scanning Calorimeter (DSC) and constant temperature ramp. Uses
   Uses Vyazkovin nonlinear method described by Sbirrazzuoli, 2009.

   Code by Craig Lawrie, (c) 2018.

   Reference:
   Sbirrazzuoli, Nicolas, et al. "Integral, differential and advanced isoconversional methods: complex
   mechanisms and isothermal predicted conversion–time curves." Chemometrics and Intelligent Laboratory
   Systems 96.2 (2009): 219-226.
"""

import sys
import numpy as np     # installed with matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from scipy.integrate import cumtrapz, trapz
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.legend_handler import HandlerLine2D

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
    #plt.figure(1)
    #plt.plot(ds1['t'], ds1['Value'])
    #plt.plot(ds2['t'], ds2['Value'])
    #plt.plot(ds3['t'], ds3['Value'])
    #plt.plot(ds4['t'], ds4['Value'])
    #plt.plot(ds5['t'], ds5['Value'])
    #plt.show()

    # Calculate total heat of reaction by integrating heat flux over time

    heat1 = cumtrapz(ds1['Value'],ds1['t'],initial = 0)
    heat2 = cumtrapz(ds2['Value'],ds2['t'],initial = 0)
    heat3 = cumtrapz(ds3['Value'],ds3['t'],initial = 0)
    heat4 = cumtrapz(ds4['Value'],ds4['t'],initial = 0)
    heat5 = cumtrapz(ds5['Value'],ds5['t'],initial = 0)

    # Plot total heat of reaction
    #plt.figure(2)
    #plt.plot(ds1['t'], heat1)
    #plt.plot(ds2['t'], heat2)
    #plt.plot(ds3['t'], heat3)
    #plt.plot(ds4['t'], heat4)
    #plt.plot(ds5['t'], heat5)
    #plt.show()

    # Calculate degree of cure
    degree_of_cure1 = heat1 / heat1[-1]
    degree_of_cure2 = heat2 / heat2[-1]
    degree_of_cure3 = heat3 / heat3[-1]
    degree_of_cure4 = heat4 / heat4[-1]
    degree_of_cure5 = heat5 / heat5[-1]

    ## Plot degree of cure
    #plt.figure(3)
    #plt.plot(ds1['t'], degree_of_cure1)
    #plt.plot(ds2['t'], degree_of_cure2)
    #plt.plot(ds3['t'], degree_of_cure3)
    #plt.plot(ds4['t'], degree_of_cure4)
    #plt.plot(ds5['t'], degree_of_cure5)
    #plt.show()

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

    # Plot rate of cure vs degree of cure AND temperature
    fig = plt.figure(4)
    ax = fig.add_subplot(111,projection='3d')
    ax.scatter(ds1['Ts'] * 9 / 5 + 32, degree_of_cure1 * 100, rate_of_cure1 * 100)
    ax.scatter(ds2['Ts'] * 9 / 5 + 32, degree_of_cure2 * 100, rate_of_cure2 * 100)
    ax.scatter(ds3['Ts'] * 9 / 5 + 32, degree_of_cure3 * 100, rate_of_cure3 * 100)
    ax.scatter(ds4['Ts'] * 9 / 5 + 32, degree_of_cure4 * 100, rate_of_cure4 * 100)
    ax.scatter(ds5['Ts'] * 9 / 5 + 32, degree_of_cure5 * 100, rate_of_cure5 * 100)
    ax.set_xlabel('Temperature (deg F)')
    ax.set_ylabel('Conversion (%)')
    ax.set_zlabel('Rate of Conversion (%/s)')
    ax.set_title('Validation Input: Rate of Cure From Dataset')
    plt.show()

    # Minimize the cost function Phi = sum(sum(J(i)/J(j))) for i != j
    # J(i) = \int_t_{_{\alpha-\Delta\alpha}}^{t_\alpha}\exp(-E/{R T_i(t))} dt
    # by varying E.
    # Whole cost function in LaTeX:
    #    \Phi(E_\alpha) = \sum_i\sum_{j\ne
    #    i}\frac{{\int_t_{_{\alpha-\Delta\alpha}}
    #      ^{t_\alpha}\exp(-E_\alpha/{R T_i(t_\alpha))}
    #      dt}}{\int_t_{_{\alpha-\Delta\alpha}}
    #      ^{t_\alpha}\exp(-E_\alpha/{R T_j(t_\alpha))} dt}
    # Repeat for a reasonable sampling of alpha.
    phi = 0.0
    J_i = 0.0
    J_j = 0.0
    runs = 5
    R = 8.314 # J/ mol.  K.
    d_alpha = 0.01
    alpha_range = np.arange(d_alpha,1.0,d_alpha)
    E_alpha = np.zeros(alpha_range.shape)

    # Set time and temperature range local to alpha value for integral
    temp_local = np.zeros([5,2])
    time_local = np.zeros([5,2])

    for alpha_step in range(alpha_range.size):
        # Reset minimization parameters
        phi_min = sys.float_info.max
        E_phi_min = 0
        E_min = 0
        E_max = 200000
        E_steps = 50

        # Set local time and temperature range for alpha
        time_local[0,0] = np.interp(np.max((alpha_range[alpha_step] - d_alpha, 0.001)), degree_of_cure1, ds1['t'])
        temp_local[0,0] = np.interp(np.max((alpha_range[alpha_step] - d_alpha, 0.001)), degree_of_cure1, ds1['Ts']) + 273.15
        time_local[1,0] = np.interp(np.max((alpha_range[alpha_step] - d_alpha, 0.001)), degree_of_cure2, ds2['t'])
        temp_local[1,0] = np.interp(np.max((alpha_range[alpha_step] - d_alpha, 0.001)), degree_of_cure2, ds2['Ts']) + 273.15
        time_local[2,0] = np.interp(np.max((alpha_range[alpha_step] - d_alpha, 0.001)), degree_of_cure3, ds3['t'])
        temp_local[2,0] = np.interp(np.max((alpha_range[alpha_step] - d_alpha, 0.001)), degree_of_cure3, ds3['Ts']) + 273.15
        time_local[3,0] = np.interp(np.max((alpha_range[alpha_step] - d_alpha, 0.001)), degree_of_cure4, ds4['t'])
        temp_local[3,0] = np.interp(np.max((alpha_range[alpha_step] - d_alpha, 0.001)), degree_of_cure4, ds4['Ts']) + 273.15
        time_local[4,0] = np.interp(np.max((alpha_range[alpha_step] - d_alpha, 0.001)), degree_of_cure5, ds5['t'])
        temp_local[4,0] = np.interp(np.max((alpha_range[alpha_step] - d_alpha, 0.001)), degree_of_cure5, ds5['Ts']) + 273.15
        time_local[0,1] = np.interp(alpha_range[alpha_step], degree_of_cure1, ds1['t'])
        temp_local[0,1] = np.interp(alpha_range[alpha_step], degree_of_cure1, ds1['Ts']) + 273.15
        time_local[1,1] = np.interp(alpha_range[alpha_step], degree_of_cure2, ds2['t'])
        temp_local[1,1] = np.interp(alpha_range[alpha_step], degree_of_cure2, ds2['Ts']) + 273.15
        time_local[2,1] = np.interp(alpha_range[alpha_step], degree_of_cure3, ds3['t'])
        temp_local[2,1] = np.interp(alpha_range[alpha_step], degree_of_cure3, ds3['Ts']) + 273.15
        time_local[3,1] = np.interp(alpha_range[alpha_step], degree_of_cure4, ds4['t'])
        temp_local[3,1] = np.interp(alpha_range[alpha_step], degree_of_cure4, ds4['Ts']) + 273.15
        time_local[4,1] = np.interp(alpha_range[alpha_step], degree_of_cure5, ds5['t'])
        temp_local[4,1] = np.interp(alpha_range[alpha_step], degree_of_cure5, ds5['Ts']) + 273.15
        for k in range(runs):
            step_size = (E_max - E_min) / E_steps
            for E_value in np.arange(E_min, E_max, step_size):
                phi = 0
                for i in range(5):
                    for j in range(5):
                        if i != j:
                            temp_knockdown = np.max((temp_local[i,0], temp_local[i,1], temp_local[j,0], temp_local[j,1])) * 2
                            J_i = 1 / 2 * (np.exp(-E_value / (R * (temp_local[i,0] - temp_knockdown))) + np.exp(-E_value / (R * (temp_local[i,1] - temp_knockdown)))) \
                                              * (time_local[i,1] - time_local[i,0])
                            J_j = 1 / 2 * (np.exp(-E_value / (R * (temp_local[j,0] - temp_knockdown))) + np.exp(-E_value / (R * (temp_local[j,1] - temp_knockdown)))) \
                                              * (time_local[j,1] - time_local[j,0])
                            
                            phi += J_i / J_j
                if phi < phi_min:
                    phi_min = phi
                    E_phi_min = E_value
            E_min = E_phi_min - step_size
            E_max = E_phi_min + step_size
        E_alpha[alpha_step] = E_phi_min
        
    

    # Plot Activation Energy
    #alpha_range = alpha_range[1:]
    #E_alpha = E_alpha[1:]
    plt.figure(5)
    plt.plot(alpha_range * 100, E_alpha / 1000)
    plt.title('Isoconversional Activation Energy for Sample Material')
    plt.xlabel('Conversion (%)')
    plt.ylabel('Activation Energy (kJ/mol)')
    plt.xlim([0,100])
    plt.show()

    # Calculate the time required to progress to a certain degree of cure at a
    # particular isothermal temperature
    # t_\alpha = [\beta \exp(-E_\alpha / {RT_{iso}})]^{-1} \int_0^{T_\alpha}
    # \exp(-E_\alpha / {RT}) dT
    temp_iso_degF_min = 200
    temp_iso_degF_max = 500
    num_temp_steps = 9

    temps_iso_degF = np.linspace(temp_iso_degF_min,temp_iso_degF_max,num_temp_steps)
    wframeX, wframeY = np.meshgrid(temps_iso_degF,alpha_range*100)
    wframeZ = np.zeros(wframeY.shape)

    plt.figure(6)

    for temp_step in range(num_temp_steps):

        temp_iso_degF = temps_iso_degF[temp_step]
        temp_iso_K = (temp_iso_degF - 32) * 5 / 9 + 273.15

        time_axis = np.zeros(alpha_range.shape)
        time_axis2 = np.zeros(alpha_range.shape)
        time_axis3 = np.zeros(alpha_range.shape)
        time_axis4 = np.zeros(alpha_range.shape)
        time_axis5 = np.zeros(alpha_range.shape)

        beta = ((ds1.Ts[2] - ds1.Ts[1]) / (ds1.t[2] - ds1.t[1]))
        for alpha_step in range(1,alpha_range.size):
            temp_min = np.interp(alpha_range[alpha_step - 1],degree_of_cure1,ds1['Ts']) + 273.15
            temp_max = np.interp(alpha_range[alpha_step],degree_of_cure1,ds1['Ts']) + 273.15
            time_axis[alpha_step] = 1 / (2 * beta) * (np.exp(-E_alpha[alpha_step] / (R * temp_max) + E_alpha[alpha_step] / (R * temp_iso_K)) + np.exp(-E_alpha[alpha_step - 1] / (R * temp_min) + E_alpha[alpha_step] / (R * temp_iso_K))) * (temp_max - temp_min)
        time_axis = np.cumsum(time_axis)

        beta = ((ds2.Ts[2] - ds2.Ts[1]) / (ds2.t[2] - ds2.t[1]))
        for alpha_step in range(1,alpha_range.size):
            temp_min = np.interp(alpha_range[alpha_step - 1],degree_of_cure2,ds2['Ts']) + 273.15
            temp_max = np.interp(alpha_range[alpha_step],degree_of_cure2,ds2['Ts']) + 273.15
            time_axis2[alpha_step] = 1 / (2 * beta) * (np.exp(-E_alpha[alpha_step] / (R * temp_max) + E_alpha[alpha_step] / (R * temp_iso_K)) + np.exp(-E_alpha[alpha_step - 1] / (R * temp_min) + E_alpha[alpha_step] / (R * temp_iso_K))) * (temp_max - temp_min)
        time_axis2 = np.cumsum(time_axis2)

        beta = ((ds3.Ts[2] - ds3.Ts[1]) / (ds3.t[2] - ds3.t[1]))
        for alpha_step in range(1,alpha_range.size):
            temp_min = np.interp(alpha_range[alpha_step - 1],degree_of_cure3,ds3['Ts']) + 273.15
            temp_max = np.interp(alpha_range[alpha_step],degree_of_cure3,ds3['Ts']) + 273.15
            time_axis3[alpha_step] = 1 / (2 * beta) * (np.exp(-E_alpha[alpha_step] / (R * temp_max) + E_alpha[alpha_step] / (R * temp_iso_K)) + np.exp(-E_alpha[alpha_step - 1] / (R * temp_min) + E_alpha[alpha_step] / (R * temp_iso_K))) * (temp_max - temp_min)
        time_axis3 = np.cumsum(time_axis3)

        beta = ((ds4.Ts[2] - ds4.Ts[1]) / (ds4.t[2] - ds4.t[1]))
        for alpha_step in range(1,alpha_range.size):
            temp_min = np.interp(alpha_range[alpha_step - 1],degree_of_cure4,ds4['Ts']) + 273.15
            temp_max = np.interp(alpha_range[alpha_step],degree_of_cure4,ds4['Ts']) + 273.15
            time_axis4[alpha_step] = 1 / (2 * beta) * (np.exp(-E_alpha[alpha_step] / (R * temp_max) + E_alpha[alpha_step] / (R * temp_iso_K)) + np.exp(-E_alpha[alpha_step - 1] / (R * temp_min) + E_alpha[alpha_step] / (R * temp_iso_K))) * (temp_max - temp_min)
        time_axis4 = np.cumsum(time_axis4)

        beta = ((ds5.Ts[2] - ds5.Ts[1]) / (ds5.t[2] - ds5.t[1]))
        for alpha_step in range(1,alpha_range.size):
            temp_min = np.interp(alpha_range[alpha_step - 1],degree_of_cure5,ds5['Ts']) + 273.15
            temp_max = np.interp(alpha_range[alpha_step],degree_of_cure5,ds5['Ts']) + 273.15
            time_axis5[alpha_step] = 1 / (2 * beta) * (np.exp(-E_alpha[alpha_step] / (R * temp_max) + E_alpha[alpha_step] / (R * temp_iso_K)) + np.exp(-E_alpha[alpha_step - 1] / (R * temp_min) + E_alpha[alpha_step] / (R * temp_iso_K))) * (temp_max - temp_min)
        time_axis5 = np.cumsum(time_axis5)

        time_axis = (time_axis + time_axis2 + time_axis3 + time_axis4 + time_axis5) / 5
        line = plt.plot(time_axis, alpha_range * 100, label = 'Temp ' + str(round(temp_iso_degF)) + ' deg F')
        plt.legend()

        for alpha_step in range(1,alpha_range.size):
            wframeZ[alpha_step, temp_step] = (alpha_range[alpha_step]-alpha_range[alpha_step-1])/(time_axis[alpha_step]-time_axis[alpha_step-1])*100
            if wframeZ[alpha_step, temp_step]>2: wframeZ[alpha_step, temp_step]=float('NaN')


    # Plot Time at Temp (still on fig 6)
    plt.title('Predicted Isothermal Cure Times')
    plt.xlabel('Time (s)')
    plt.ylabel('Conversion (%)')
    plt.xlim([0,3600])
    plt.show()
    
    plt.ioff()

    # Plot validation against cure rate from input data
    fig = plt.figure(7)
    ax = fig.add_subplot(111,projection='3d')
    ax.plot_wireframe(wframeX,wframeY,wframeZ)
    ax.scatter(ds1['Ts'] * 9 / 5 + 32, degree_of_cure1 * 100, rate_of_cure1 * 100)
    ax.scatter(ds2['Ts'] * 9 / 5 + 32, degree_of_cure2 * 100, rate_of_cure2 * 100)
    ax.scatter(ds3['Ts'] * 9 / 5 + 32, degree_of_cure3 * 100, rate_of_cure3 * 100)
    ax.scatter(ds4['Ts'] * 9 / 5 + 32, degree_of_cure4 * 100, rate_of_cure4 * 100)
    ax.scatter(ds5['Ts'] * 9 / 5 + 32, degree_of_cure5 * 100, rate_of_cure5 * 100)
    ax.set_zlim([0,2])
    ax.set_xlabel('Temperature (deg F)')
    ax.set_ylabel('Conversion (%)')
    ax.set_zlabel('Rate of Conversion (%/s)')
    ax.set_title('Validation Output: Rate of Cure From Isothermal Predictions')
    plt.show()



main()