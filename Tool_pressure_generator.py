import sys
sys.path.append('/data1/caldas/Pytmosph3R/PyRouts/')
from pyfunction import ratio
from pyfunction import interp2olation
from pyfunction import interpolation
from pyfunction import interpolation_multi
import scipy.integrate as integrate
from pyconstant import *
import numpy as np
import matplotlib.pyplot as plt
import pickle
from scipy.integrate import quad

def pressure_gravity(T_iso,Rp,M,g0,P_surf,P_h,n_layers,n_species,number,x_ratio_species,Integration,Layers,Plotting) :

    lim_alt = R_gp*T_iso/(M*g0)*np.log(P_surf/P_h)
    rupt_alt, h = 0.e+0, lim_alt
    delta_z = h/np.float(n_layers)
    r_step, x_step, theta_number, reso_alt = delta_z, delta_z, 1, int(h/1000)
    z_lim = int(lim_alt/delta_z)
    z_reso = int(h/delta_z) + 1
    z_array = np.arange(h/np.float(delta_z)+1)*float(delta_z)
    print 'pressure_gravity says : Scale height : ',R_gp*T_iso/(M*g0),' m'
    data = np.zeros((number,n_layers+2))
    data[0,0], data[1,0], data[2:2+n_species.size,0], data[number-1,0] = P_surf, T_iso, x_ratio_species, M
    data_cc = np.zeros((number,n_layers+2))
    data_cc[0,0], data_cc[1,0], data_cc[2:2+n_species.size,0], data_cc[number-1,0] = P_surf, T_iso, x_ratio_species, M
    alt_array = np.zeros(n_layers+2)
    alt_array[0] = 0

    for i_l in range(1,n_layers+2) :
        if i_l != n_layers+1 :
            z = (i_l - 0.5)*delta_z
        else :
            z = (i_l - 1)*delta_z
        alt_array[i_l] = z
        if i_l == 1 or i_l == n_layers+1 :
            if Integration == True :
                data[0,i_l] = P_surf*np.exp(-M*g0/(R_gp*T_iso)*z)
            if Layers == True :
                data_cc[0,i_l] = data_cc[0,i_l-1]*np.exp(-M*g0/(R_gp*T_iso)*delta_z/2.)
        else :
            if Integration == True :
                data[0,i_l] = P_surf*np.exp(-M*g0/(R_gp*T_iso)*z)
            if Layers == True :
                data_cc[0,i_l] = data_cc[0,i_l-1]*np.exp(-M*g0/(R_gp*T_iso)*delta_z)

        data[1,i_l], data[2:2+n_species.size,i_l], data[number-1,i_l] = T_iso, x_ratio_species, M
        data_cc[1,i_l], data_cc[2:2+n_species.size,i_l], data_cc[number-1,i_l] = T_iso, x_ratio_species, M

    if Plotting == True :
        if Integration == True :
            plt.plot(alt_array,data[0],'+-g',linewidth=4, ms=12, label='Gravity constant profile (integration)')
        if Layers == True :
            plt.plot(alt_array,data_cc[0], '+-r',linewidth=2, ms=6, label='Gravity constant profile')

    return data, data_cc, alt_array, h, delta_z, r_step, x_step, theta_number, reso_alt, z_lim, z_reso, z_array

########################################################################################################################

def pressure_scaling(T_iso,Rp,M,g0,P_surf,P_h,n_layers,n_species,number,x_ratio_species,Integration,Layers,Plotting) :

    p = np.zeros(n_layers+2)
    p[0] = np.log(P_surf)
    print 'pressure_scaling says : Scale height : ',R_gp*T_iso/(M*g0),' m'
    data = np.zeros((number,n_layers+2))
    data[0,0], data[1,0], data[2:2+n_species.size,0], data[number-1,0] = P_surf, T_iso, x_ratio_species, M
    alt_array = np.zeros(n_layers+2)
    alt_array[0] = 0
    alt_array_cc = np.zeros(n_layers+2)
    alt_array_cc[0] = 0
    g_array = np.zeros(n_layers+2)
    g_array[0] = g0
    H_array = np.zeros(n_layers+2)
    H_array[0] = k_B*T_iso/(M*AMU*g0)*1.e-3

    for i_p in range(1,n_layers+2) :
        if i_p != n_layers + 1 :
            p[i_p] = np.log(P_surf) - (i_p - 0.5)*np.log(P_surf/P_h)/np.float(n_layers+1)
        else :
            p[i_p] = np.log(P_h)
        data[0,i_p] = np.exp(p[i_p])
        if Integration == True :
            alt_array[i_p] = (np.log(P_surf)-p[i_p])/(M*g0/(R_gp*T_iso)-1./Rp*(np.log(P_surf)-p[i_p]))
            #print 'integration', alt_array[i_p]
            g_array[i_p] = g0*1./(1+alt_array[i_p]/Rp)**2
            H_array[i_p] = k_B*T_iso/(M*AMU*g_array[i_p])*1.e-3
        if Layers == True :
            alt_array_cc[i_p] = alt_array_cc[i_p-1] + (p[i_p-1]-p[i_p])*\
                (1+alt_array_cc[i_p-1]/Rp)/(M*g0/(R_gp*T_iso)*1./(1+alt_array_cc[i_p-1]/Rp)-1./Rp*(p[i_p-1]-p[i_p]))
            #print 'line by line', alt_array_cc[i_p]
            g_array[i_p] = g0*1./(1+alt_array_cc[i_p]/Rp)**2
            H_array[i_p] = k_B*T_iso/(M*AMU*g_array[i_p])*1.e-3
        if i_p == n_layers+1 :
            if Integration == True :
                lim_alt = alt_array[i_p]
            if Layers == True :
                lim_alt = alt_array_cc[i_p]

        data[1,i_p], data[2:2+n_species.size,i_p], data[number-1,i_p] = T_iso, x_ratio_species, M

    rupt_alt, h = 0.e+0, lim_alt
    delta_z = h/np.float(n_layers)
    r_step, x_step, theta_number, reso_alt = delta_z, delta_z, 1, int(h/1000)
    z_lim = int(lim_alt/delta_z)
    z_reso = int(h/delta_z) + 1
    z_array = np.arange(h/np.float(delta_z)+1)*float(delta_z)
    if Plotting == True :
        if Integration == True :
            plt.plot(alt_array,data[0],'+--r',linewidth=4, ms=16, label='Hydrostatic profile (integration)')
        if Layers == True :
            plt.plot(alt_array_cc,data[0], '+--b', linewidth=4, ms=16, label='Hydrostatic profile')

    return data, alt_array, alt_array_cc, H_array, g_array, h, delta_z, r_step, x_step, theta_number, reso_alt, z_lim, z_reso, z_array

########################################################################################################################

def pressure_scaling_taurex(T_iso,Rp,M,Mp,P_surf,P_h,n_layers,n_species,number,x_ratio_species,Plotting,Script=True) :

    press_exp = np.linspace(np.log(P_h), np.log(P_surf), n_layers+1)
    pressure_levels = np.exp(press_exp)[::-1]
    data = np.zeros((number,n_layers+2))
    data[0,0], data[1,0], data[2:2+n_species.size,0], data[number-1,0] = P_surf, T_iso, x_ratio_species, M
    data[0,1:n_layers+1] = np.power(10,np.log10(pressure_levels[:-1])+\
                    np.diff(np.log10(pressure_levels))/2.)
    data[0,n_layers+1] = P_h
    alt_array = np.zeros(n_layers+1)
    alt_array[0] = 0
    g_array = np.zeros(n_layers+1)
    g_array[0] = G*Mp/(Rp**2)
    H_array = np.zeros(n_layers+1)
    H_array[0] = k_B*T_iso/(M*AMU*g_array[0])*1.e-3
    if Script == True :
        print 'pressure_scaling_taurex says : Scale height : ',R_gp*T_iso/(M*g_array[0]),' m'

    for i_l in range(1,n_layers) :

        delta_r = (-1.)*H_array[i_l-1]*np.log(pressure_levels[i_l]/pressure_levels[i_l-1])
        alt_array[i_l] = alt_array[i_l-1] + delta_r
        g_array[i_l] = G*Mp/(Rp+alt_array[i_l])**2
        H_array[i_l] = k_B*T_iso/(M*AMU*g_array[i_l])*1.e-3

        data[1,i_l], data[2:2+n_species.size,i_l], data[number-1,i_l] = T_iso, x_ratio_species, M

    alt_array[n_layers] = alt_array[n_layers-1] + delta_r
    g_array[n_layers] = G*Mp/(Rp+alt_array[n_layers])**2
    H_array[n_layers] = k_B*T_iso/(M*AMU*g_array[n_layers])*1.e-3
    data[1,n_layers], data[2:2+n_species.size,n_layers], data[number-1,n_layers] = T_iso, x_ratio_species, M
    data[1,n_layers+1], data[2:2+n_species.size,n_layers+1], data[number-1,n_layers+1] = T_iso, x_ratio_species, M
    lim_alt = alt_array[n_layers]
    rupt_alt, h = 0.e+0, lim_alt
    delta_z = h/np.float(n_layers)
    r_step, x_step, theta_number, reso_alt = delta_z, delta_z, 1, int(h/1000)
    z_lim = int(lim_alt/delta_z)
    z_reso = int(h/delta_z) + 1
    z_array = np.arange(h/np.float(delta_z)+1)*float(delta_z)
    alt_mid = np.zeros(n_layers)
    for i_w in range(n_layers) :
        if i_w != n_layers-1 :
            alt_mid[i_w] = (alt_array[i_w+1] + alt_array[i_w])/2.
        else :
            alt_mid[i_w] = alt_array[i_w] + (alt_array[i_w] - alt_array[i_w-1])/2.

    if Plotting == True :
        plt.plot(alt_mid,data[0,1:n_layers+1],'+--g',linewidth=4, ms=16, label='Taurex writing')

    return data, alt_array, alt_mid, H_array, g_array, h, delta_z, r_step, x_step, theta_number, reso_alt, z_lim, z_reso, z_array


########################################################################################################################

def pressure_gravity_evolution(T_iso,Rp,M,g0,P_surf,P_h,n_layers,n_species,number,x_ratio_species,Integration,Layers,Other,Taurdu,Plotting) :

    lim_alt = R_gp*T_iso/(M*g0)*np.log(P_surf/P_h)*1./(1. + R_gp*T_iso/(Rp*M*g0)*np.log(P_h/P_surf))
    print lim_alt
    rupt_alt, h = 0.e+0, lim_alt
    delta_z = h/np.float(n_layers)
    r_step, x_step, theta_number, reso_alt = delta_z, delta_z, 1, int(h/1000)
    z_lim = int(lim_alt/delta_z)
    z_reso = int(h/delta_z) + 1
    z_array = np.arange(h/np.float(delta_z)+1)*float(delta_z)
    print 'pressure_gravity_evolution says : Scale height : ',R_gp*T_iso/(M*g0),' m'
    data = np.zeros((number,n_layers+2))
    data_cc = np.zeros((number,n_layers+2))
    data_cc_o = np.zeros((number,n_layers+2))
    data_cc_t = np.zeros((number,n_layers+2))
    data[0,0], data[1,0], data[2:2+n_species.size,0], data[number-1,0] = P_surf, T_iso, x_ratio_species, M
    data_cc[0,0], data_cc[1,0], data_cc[2:2+n_species.size,0], data_cc[number-1,0] = P_surf, T_iso, x_ratio_species, M
    data_cc_o[0,0], data_cc_o[1,0], data_cc_o[2:2+n_species.size,0], data_cc_o[number-1,0] = P_surf, T_iso, x_ratio_species, M
    data_cc_t[0,0], data_cc_t[1,0], data_cc_t[2:2+n_species.size,0], data_cc_t[number-1,0] = P_surf, T_iso, x_ratio_species, M
    g_array = np.zeros(n_layers+2)
    alt_array = np.zeros(n_layers+2)
    H_array = np.zeros(n_layers+2)
    alt_array[0] = 0
    g_array[0] = g0
    H_array[0] = k_B*T_iso/(M*AMU*g0)*1.e-3

    for i_l in range(1,n_layers+2) :
        if i_l != n_layers+1 :
            z = (i_l - 0.5)*delta_z
        else :
            z = (i_l - 1)*delta_z
        alt_array[i_l] = z
        g_array[i_l] = g0*1./(1+z/Rp)**2
        H_array[i_l] = k_B*T_iso/(M*AMU*g_array[i_l])*1.e-3
        if i_l == 1 or i_l == n_layers+1 :
            if Integration == True :
                data[0,i_l] = P_surf*np.exp(-M*g0/(R_gp*T_iso)*z/(1.+z/Rp))
            if Layers == True :
                data_cc[0,i_l] = data_cc[0,i_l-1]*np.exp(-M*g0/(R_gp*T_iso)*delta_z/2.*1./((1+(z-0.5*delta_z)/Rp)**2\
                                                                   *(1+0.5*delta_z/(Rp+z-0.5*delta_z))))
            if Other == True :
                data_cc_o[0,i_l] = data_cc_o[0,i_l-1]*np.exp(-M*g0/(R_gp*T_iso)*0.5*delta_z*1./((1+(z-0.5*delta_z)/Rp)*(1+z/(Rp))))
            if Taurdu == True :
                data_cc_t[0,i_l] = data_cc_t[0,i_l-1]*np.exp(-M*g_array[i_l-1]/(R_gp*T_iso)*0.5*delta_z)
        else :
            if Integration == True :
                data[0,i_l] = P_surf*np.exp(-M*g0/(R_gp*T_iso)*z/(1.+z/Rp))
            if Layers == True :
                data_cc[0,i_l] = data_cc[0,i_l-1]*np.exp(-M*g0/(R_gp*T_iso)*delta_z*1./((1+(z-delta_z)/Rp)**2\
                                                                                *(1+delta_z/(Rp+z-delta_z))))
            if Other == True :
                data_cc_o[0,i_l] = data_cc_o[0,i_l-1]*np.exp(-M*g0/(R_gp*T_iso)*delta_z*1./((1+(z-delta_z)/Rp)*(1+z/(Rp))))
            if Taurdu == True :
                data_cc_t[0,i_l] = data_cc_t[0,i_l-1]*np.exp(-M*g_array[i_l-1]/(R_gp*T_iso)*delta_z)

        data[1,i_l], data[2:2+n_species.size,i_l], data[number-1,i_l] = T_iso, x_ratio_species, M
        data_cc[1,i_l], data_cc[2:2+n_species.size,i_l], data_cc[number-1,i_l] = T_iso, x_ratio_species, M
        data_cc_o[1,i_l], data_cc_o[2:2+n_species.size,i_l], data_cc_o[number-1,i_l] = T_iso, x_ratio_species, M
        data_cc_t[1,i_l], data_cc_t[2:2+n_species.size,i_l], data_cc_t[number-1,i_l] = T_iso, x_ratio_species, M

    if Plotting == True :
        if Taurdu == True :
            plt.plot(alt_array,data_cc_t[0], '+-y', linewidth=6, ms=24, label='Hydrostatic level profil')
        if Other == True :
            plt.plot(alt_array,data_cc_o[0], '+-c', linewidth=6, ms=24)
        if Integration == True :
            plt.plot(alt_array,data[0],'+-b',linewidth=4, ms=16, label='Hydrostatic profile (integration)')
        if Layers == True :
            plt.plot(alt_array,data_cc[0], '+-r', linewidth=2, ms=8, label='Hydrostatic profile')
            #np.save('/Users/caldas/Desktop/Pytmosph3R/Tools/Output_Tools/data.npy',data_cc)
            #np.save('/Users/caldas/Desktop/Pytmosph3R/Tools/Output_Tools/alt.npy',alt_array)

    return data, data_cc, data_cc_o, data_cc_t, alt_array, H_array, g_array, h, delta_z, r_step, x_step, theta_number, reso_alt, z_lim, z_reso, z_array

########################################################################################################################

def pressure_evolution_mass(T_iso,Rp,M,g0,P_surf,P_h,n_layers,n_species,number,x_ratio_species,Integration,Layers,Plotting) :

    alp_alt = R_gp*T_iso/(M*g0)*np.log(P_h/P_surf)
    lim_alt = -alp_alt/(1. + alp_alt/Rp)
    rupt_alt, h = 0.e+0, lim_alt
    delta_z = h/np.float(n_layers)
    r_step, x_step, theta_number, reso_alt = delta_z, delta_z, 1, int(h/1000)
    z_lim = int(lim_alt/delta_z)
    z_reso = int(h/delta_z) + 1
    z_array = np.arange(h/np.float(delta_z)+1)*float(delta_z)
    print 'pressure_evolution_mass says : Scale height : ',R_gp*T_iso/(M*g0),' m'
    data = np.zeros((number,n_layers+2))
    data_cc = np.zeros((number,n_layers+2))
    data[0,0], data[1,0], data[2:2+n_species.size,0], data[number-1,0] = P_surf, T_iso, x_ratio_species, M
    data_cc[0,0], data_cc[1,0], data_cc[2:2+n_species.size,0], data_cc[number-1,0] = P_surf, T_iso, x_ratio_species, M
    excess = 0
    g_plus = 0
    alt_array = np.zeros(n_layers+2)
    g_array = np.zeros(n_layers+2)
    alt_array[0] = 0
    g_array[0] = g0
    H_array = np.zeros(n_layers+2)
    H_array[0] = k_B*T_iso/(M*AMU*g0)*1.e-3

    for i_l in range(1,n_layers+2) :
        if i_l != n_layers+1 :
            z = (i_l - 0.5)*delta_z
        else :
            z = (i_l - 1)*delta_z
        alt_array[i_l] = z
        g_array[i_l] = g0*1./(1+z/Rp)**2 + g_plus
        H_array[i_l] = k_B*T_iso/(M*AMU*g_array[i_l])*1.e-3
        if i_l == 1 or i_l == n_layers+1 :
            if i_l == 1 :
                if Integration == True :
                    Int = integrate.quad(lambda x:-M/(R_gp*T_iso)*(g0*(1./(1+x/Rp)**2)+\
                                G*M/(R_gp*T_iso)*4/3.*np.pi*P_surf*np.exp(-M*g0/(R_gp*T_iso)*z/(1+x/Rp))*((Rp+x)**3-Rp**3)/(Rp+x)**2),0,z)
                    data[0,i_l] = P_surf*np.exp(Int[0])
                if Layers == True :
                    Int = integrate.quad(lambda x:-M/(R_gp*T_iso)*((g0*(1./(1+(z-0.5*delta_z)/Rp)**2))*\
                                1./(1+(x-z+0.5*delta_z)/(Rp+z-0.5*delta_z))**2+G*M/(R_gp*T_iso)*4/3.*np.pi*data[0,i_l-1]*\
                                np.exp(-M*g0/(R_gp*T_iso)*x/(1+x/Rp))*((Rp+x)**3-(Rp+z-0.5*delta_z)**3)/(Rp+x)**2),z-0.5*delta_z,z)
                    data_cc[0,i_l] = data_cc[0,i_l-1]*np.exp(Int[0])
                    excess = data_cc[0,i_l-1]*M/(R_gp*T_iso)*4/3*np.pi*((Rp+z)**3-(Rp+z-0.5*delta_z)**3)
            else :
                g_plus += excess*G/(Rp+z)**2
                if Integration == True :
                    Int = integrate.quad(lambda x:-M/(R_gp*T_iso)*(g0*(1./(1+x/Rp)**2)+\
                                G*M/(R_gp*T_iso)*4/3.*np.pi*P_surf*np.exp(-M*g0/(R_gp*T_iso)*z/(1+x/Rp))*((Rp+x)**3-Rp**3)/(Rp+x)**2),0,z)
                    data[0,i_l] = P_surf*np.exp(Int[0])
                if Layers == True :
                    Int = integrate.quad(lambda x:-M/(R_gp*T_iso)*((g0*(1./(1+(z-0.5*delta_z)/Rp)**2)+g_plus)*\
                            1./(1+(x-z+0.5*delta_z)/(Rp+z-0.5*delta_z))**2+G*M/(R_gp*T_iso)*4/3.*np.pi*data[0,i_l-1]*\
                            np.exp(-M*g0/(R_gp*T_iso)*x/(1+x/Rp))*((Rp+x)**3-(Rp+z-0.5*delta_z)**3)/(Rp+x)**2),z-0.5*delta_z,z)
                    data_cc[0,i_l] = data_cc[0,i_l-1]*np.exp(Int[0])
        else :
            if Integration == True :
                Int = integrate.quad(lambda x:-M/(R_gp*T_iso)*(g0*(1./(1+x/Rp)**2)+\
                            G*M/(R_gp*T_iso)*4/3.*np.pi*P_surf*np.exp(-M*g0/(R_gp*T_iso)*\
                            z/(1+x/Rp))*((Rp+x)**3-Rp**3)/(Rp+x)**2),0,z)
                data[0,i_l] = P_surf*np.exp(Int[0])
            if Layers == True :
                g_plus += excess*G/(Rp+z-0.5*delta_z)**2
                Int = integrate.quad(lambda x:-M/(R_gp*T_iso)*((g0*(1./(1+(z-delta_z)/Rp)**2)+g_plus)*\
                            1./(1+(x-z+delta_z)/(Rp+z-delta_z))**2+G*M/(R_gp*T_iso)*4/3.*np.pi*data[0,i_l-1]*\
                            np.exp(-M*g0/(R_gp*T_iso)*x/(1+x/Rp))*((Rp+x)**3-(Rp+z-delta_z)**3)/(Rp+x)**2),z-delta_z,z)
                data_cc[0,i_l] = data_cc[0,i_l-1]*np.exp(Int[0])
                excess = data_cc[0,i_l-1]*M/(R_gp*T_iso)*4/3*np.pi*((Rp+z)**3-(Rp+z-delta_z)**3)

        data[1,i_l], data[2:2+n_species.size,i_l], data[number-1,i_l] = T_iso, x_ratio_species, M
        data_cc[1,i_l], data_cc[2:2+n_species.size,i_l], data_cc[number-1,i_l] = T_iso, x_ratio_species, M
    if Plotting == True :
        if Integration == True :
            plt.plot(alt_array,data[0],'+-m', linewidth=4, ms=16, label='Pytmosph3R profile (integration)')
        if Layers == True :
            plt.plot(alt_array,data_cc[0], '+-y', linewidth=2, ms=8, label='Pytmosph3R profile')
            print H_array

    return data, data_cc, alt_array, h, H_array, g_array, delta_z, r_step, x_step, theta_number, reso_alt, z_lim, z_reso, z_array

########################################################################################################################

def instance_extract(taurex_file) :

    instance = pickle.load(open(taurex_file,'rb'))
    P_pro = instance['data']['altitude_profile']
    alt = P_pro[:,1]
    alt_mid = np.zeros(alt.size)
    pressure = P_pro[:,0]
    for i_w in range(alt.size) :
        if i_w != alt.size-1 :
            alt_mid[i_w] = (alt[i_w+1] + alt[i_w])/2.
        else :
            alt_mid[i_w] = alt[i_w] + (alt[i_w] - alt[i_w-1])/2.
    plt.plot(alt_mid,pressure,'+-k',linewidth=3, ms=12, label='Taurex profile')

    T = instance['data']['temperature_profile'][:,1]
    M = instance['data']['mu_profile'][:,1]
    g = instance['data']['gravity_profile'][:,1]
    H = instance['data']['scaleheight_profile'][:,1]
    x_active = instance['data']['active_mixratio_profile'][:,1]
    x_inactive = instance['data']['inactive_mixratio_profile'][:,1]

    return alt, alt_mid, pressure, T, M, g, H, x_active, x_inactive

########################################################################################################################

def pressure_scaling_noiso_taurex(T,Rp,Mp,P_surf,P_h,n_layers,n_species,number,x_ratio_species,T_comp,P_comp,Plotting) :

    press_exp = np.linspace(np.log(P_h), np.log(P_surf), n_layers+1)
    pressure_levels = np.exp(press_exp)[::-1]
    data = np.zeros((number,n_layers+2))
    data[0,0], data[1,0] = P_surf, T[1,0]
    data[0,1:n_layers+1] = np.power(10,np.log10(pressure_levels[:-1])+\
                    np.diff(np.log10(pressure_levels))/2.)
    data[0,n_layers+1] = P_h
    alt_array = np.zeros(n_layers+1)
    alt_array[0] = 0
    g_array = np.zeros(n_layers+1)
    g_array[0] = G*Mp/(Rp**2)
    H_array = np.zeros(n_layers+1)

    data[1,1:] = np.interp(data[0,1:],T[0],T[1])
    for i_n in range(n_species.size-2):
        data[4] = interp2olation(data[0],data[1],P_comp,T_comp,x_ratio_species[i_n])
    for i_P in range(data[0].size) :
        M_species, M_result, result = ratio(n_species,data[4:2+n_species.size,i_P],IsoComp=True)
        data[2,i_P] = result[0]
        data[3,i_P] = result[1]
        data[number-1,i_P] = M_result

    H_array[0] = k_B*data[1,0]/(data[number-1,0]*AMU*g_array[0])*1.e-3
    print 'pressure_scaling_taurex says : Scale height : ',R_gp*data[1,0]/(data[number-1,0]*g_array[0]),' m'

    for i_l in range(1,n_layers) :

        delta_r = (-1.)*H_array[i_l-1]*np.log(pressure_levels[i_l]/pressure_levels[i_l-1])
        alt_array[i_l] = alt_array[i_l-1] + delta_r
        g_array[i_l] = G*Mp/(Rp+alt_array[i_l])**2
        H_array[i_l] = k_B*data[1,i_l]/(data[number-1,i_l]*AMU*g_array[i_l])*1.e-3

    alt_array[n_layers] = alt_array[n_layers-1] + delta_r
    g_array[n_layers] = G*Mp/(Rp+alt_array[n_layers])**2
    H_array[n_layers] = k_B*data[1,n_layers]/(data[number-1,n_layers]*AMU*g_array[n_layers])*1.e-3

    lim_alt = alt_array[n_layers]
    rupt_alt, h = 0.e+0, lim_alt
    delta_z = h/np.float(n_layers)
    r_step, x_step, theta_number, reso_alt = delta_z, delta_z, 1, int(h/1000)
    z_lim = int(lim_alt/delta_z)
    z_reso = int(h/delta_z) + 1
    z_array = np.arange(h/np.float(delta_z)+1)*float(delta_z)
    alt_mid = np.zeros(n_layers)
    for i_w in range(n_layers) :
        if i_w != n_layers-1 :
            alt_mid[i_w] = (alt_array[i_w+1] + alt_array[i_w])/2.
        else :
            alt_mid[i_w] = alt_array[i_w] + (alt_array[i_w] - alt_array[i_w-1])/2.

    if Plotting == True :
        plt.plot(alt_mid,data[0,1:n_layers+1],'+--g',linewidth=4, ms=16, label='Taurex writing')

    return data, alt_array, alt_mid, H_array, g_array, h, delta_z, r_step, x_step, theta_number, reso_alt, z_lim, z_reso, z_array

########################################################################################################################

def pressure_scaling_noiso_cloud_taurex(T,Rp,Mp,P_surf,P_h,n_layers,n_species,number,x_ratio_species,T_comp,P_comp,
                        c_species,gen_cond,Clouds,Plotting,Atmass=False) :

    Mp_e = Mp
    press_exp = np.linspace(np.log(P_h), np.log(P_surf), n_layers+1)
    pressure_levels = np.exp(press_exp)[::-1]
    data = np.zeros((number,n_layers+2))
    data[0,0], data[1,0] = P_surf, T[1,0]
    data[0,1:n_layers+1] = np.power(10,np.log10(pressure_levels[:-1])+\
                    np.diff(np.log10(pressure_levels))/2.)
    data[0,n_layers+1] = P_h
    alt_array = np.zeros(n_layers+1)
    alt_array[0] = 0
    g_array = np.zeros(n_layers+1)
    g_array[0] = G*Mp/(Rp**2)
    H_array = np.zeros(n_layers+1)

    data[1,1:] = np.interp(data[0,1:],T[0],T[1])
    for i_n in range(n_species.size-2):
        data[4+c_species.size],c_gr,i_gr = interp2olation(data[0],data[1],P_comp,T_comp,x_ratio_species[i_n])
    for i_P in range(data[0].size) :
        M_species, M_result, result = ratio(n_species,data[4:2+n_species.size,i_P],IsoComp=True)
        data[2+c_species.size,i_P] = result[0]
        data[3+c_species.size,i_P] = result[1]
        data[number-1,i_P] = M_result
    if Clouds == True :
        for i_c in range(c_species.size) :
            data[2+i_c],c_gr,i_gr = interp2olation(data[0],data[1],P_comp,T_comp,gen_cond[i_c])*data[number-1]/M_mole[np.where(M_n_mole == c_species[i_c])[0]]

    H_array[0] = k_B*data[1,0]/(data[number-1,0]*AMU*g_array[0])*1.e-3
    print 'pressure_scaling_taurex says : Scale height : ',R_gp*data[1,0]/(data[number-1,0]*g_array[0]),' m'

    for i_l in range(1,n_layers) :

        delta_r = (-1.)*H_array[i_l-1]*np.log(pressure_levels[i_l]/pressure_levels[i_l-1])
        alt_array[i_l] = alt_array[i_l-1] + delta_r
        if Atmass == True :
            r = alt_array[i_l] - alt_array[i_l-1]
            factor = 4*np.pi*data[number-1,i_l]*data[0,i_l-1]/(R_gp*data[1,i_l])
            alp = -data[number-1,i_l]*Mp*G/(R_gp*data[1,i_l]*(Rp+alt_array[i_l])**2)
            M_atm, l = quad(lambda z:factor*np.exp(alp*1./(1+z/(Rp+alt_array[i_l])**2))*z**2,0,r)
            if M_atm > 0.001*Mp :
                decade = np.log10(M_atm/Mp)+4
                n_sub = np.int(10**(decade/2.))
                step = (np.log10(data[0,i_l-1])-np.log10(data[0,i_l]))/np.float(n_sub)
                alt = alt_array[i_l-1]
                for i_sub in range(n_sub) :
                    H_sub = k_B*data[1,i_l]/(data[number-1,i_l]*AMU*(Mp*G/(Rp+alt)**2))*1.e-3
                    P_sub_0 = 10**(np.log10(data[0,i_l-1])-(i_sub)*step)
                    P_sub = 10**(np.log10(data[0,i_l-1])-(i_sub+1)*step)
                    delta_sub = (-1.)*H_sub*np.log(P_sub/P_sub_0)
                    alt += delta_sub
                    r = delta_sub
                    factor = 4*np.pi*data[number-1,i_l]*P_sub_0/(R_gp*data[1,i_l])
                    alp = -data[number-1,i_l]*Mp*G/(R_gp*data[1,i_l]*(Rp+alt)**2)
                    M_atm, l = quad(lambda z:factor*np.exp(alp*1./(1+z/(Rp+alt)**2))*z**2,0,r)
                    Mp += M_atm
                alt_array[i_l] = alt
            else :
                Mp += M_atm

        g_array[i_l] = G*Mp/(Rp+alt_array[i_l])**2
        H_array[i_l] = k_B*data[1,i_l]/(data[number-1,i_l]*AMU*g_array[i_l])*1.e-3

    if Mp_e != Mp :
        print 'Mass correction : %.10E kg, atmospheric mass : %.10E kg'%(Mp, Mp-Mp_e)

    alt_array[n_layers] = alt_array[n_layers-1] + delta_r
    g_array[n_layers] = G*Mp/(Rp+alt_array[n_layers])**2
    H_array[n_layers] = k_B*data[1,n_layers]/(data[number-1,n_layers]*AMU*g_array[n_layers])*1.e-3

    lim_alt = alt_array[n_layers]
    rupt_alt, h = 0.e+0, lim_alt
    delta_z = h/np.float(n_layers)
    r_step, x_step, theta_number, reso_alt = delta_z, delta_z, 1, int(h/1000)
    z_lim = int(lim_alt/delta_z)
    z_reso = int(h/delta_z) + 1
    z_array = np.arange(h/np.float(delta_z)+1)*float(delta_z)
    alt_mid = np.zeros(n_layers)
    for i_w in range(n_layers) :
        if i_w != n_layers-1 :
            alt_mid[i_w] = (alt_array[i_w+1] + alt_array[i_w])/2.
        else :
            alt_mid[i_w] = alt_array[i_w] + (alt_array[i_w] - alt_array[i_w-1])/2.

    if Plotting == True :
        plt.plot(alt_mid,data[0,1:n_layers+1],'+--g',linewidth=4, ms=16, label='Taurex writing')

    return data, alt_array, alt_mid, H_array, g_array, h, delta_z, r_step, x_step, theta_number, reso_alt, z_lim, z_reso, z_array

########################################################################################################################

def pressure_scaling_noiso_cloud_taurex_gravity(T,Rp,Mp,P_surf,P_h,n_layers,n_species,number,x_ratio_species,T_comp,P_comp,
                        c_species,gen_cond,Clouds,Plotting) :

    press_exp = np.linspace(np.log(P_h), np.log(P_surf), n_layers+1)
    pressure_levels = np.exp(press_exp)[::-1]
    data = np.zeros((number,n_layers+2))
    data[0,0], data[1,0] = P_surf, T[1,0]
    data[0,1:n_layers+1] = np.power(10,np.log10(pressure_levels[:-1])+\
                    np.diff(np.log10(pressure_levels))/2.)
    data[0,n_layers+1] = P_h
    alt_array = np.zeros(n_layers+1)
    alt_array[0] = 0
    g_array = np.zeros(n_layers+1)
    g_array[0] = G*Mp/(Rp**2)
    H_array = np.zeros(n_layers+1)

    data[1,1:] = np.interp(data[0,1:],T[0],T[1])
    for i_n in range(n_species.size-2):
        data[4+c_species.size],c_gr,i_gr = interp2olation(data[0],data[1],P_comp,T_comp,x_ratio_species[i_n])
    for i_P in range(data[0].size) :
        M_species, M_result, result = ratio(n_species,data[4:2+n_species.size,i_P],IsoComp=True)
        data[2+c_species.size,i_P] = result[0]
        data[3+c_species.size,i_P] = result[1]
        data[number-1,i_P] = M_result
    if Clouds == True :
        for i_c in range(c_species.size) :
            data[2+i_c],c_gr,i_gr = interp2olation(data[0],data[1],P_comp,T_comp,gen_cond[i_c])*data[number-1]/M_mole[np.where(M_n_mole == c_species[i_c])[0]]

    H_array[0] = k_B*data[1,0]/(data[number-1,0]*AMU*g_array[0])*1.e-3
    print 'pressure_scaling_taurex says : Scale height : ',R_gp*data[1,0]/(data[number-1,0]*g_array[0]),' m'

    for i_l in range(1,n_layers) :

        delta_r = (-1.)*H_array[i_l-1]*np.log(pressure_levels[i_l]/pressure_levels[i_l-1])
        alt_array[i_l] = alt_array[i_l-1] + delta_r
        g_array[i_l] = g_array[0]
        H_array[i_l] = k_B*data[1,i_l]/(data[number-1,i_l]*AMU*g_array[i_l])*1.e-3

    alt_array[n_layers] = alt_array[n_layers-1] + delta_r
    g_array[n_layers] = g_array[0]
    H_array[n_layers] = k_B*data[1,n_layers]/(data[number-1,n_layers]*AMU*g_array[n_layers])*1.e-3

    lim_alt = alt_array[n_layers]
    rupt_alt, h = 0.e+0, lim_alt
    delta_z = h/np.float(n_layers)
    r_step, x_step, theta_number, reso_alt = delta_z, delta_z, 1, int(h/1000)
    z_lim = int(lim_alt/delta_z)
    z_reso = int(h/delta_z) + 1
    z_array = np.arange(h/np.float(delta_z)+1)*float(delta_z)
    alt_mid = np.zeros(n_layers)
    for i_w in range(n_layers) :
        if i_w != n_layers-1 :
            alt_mid[i_w] = (alt_array[i_w+1] + alt_array[i_w])/2.
        else :
            alt_mid[i_w] = alt_array[i_w] + (alt_array[i_w] - alt_array[i_w-1])/2.

    if Plotting == True :
        plt.plot(alt_mid,data[0,1:n_layers+1],'+--g',linewidth=4, ms=16, label='Taurex writing')

    return data, alt_array, alt_mid, H_array, g_array, h, delta_z, r_step, x_step, theta_number, reso_alt, z_lim, z_reso, z_array

########################################################################################################################

def tropo_reconstruction(R_e,Mp,P_t,T_t,z_t,tropo,M,gamma,n_layers,number,fact,n_species,x_ratio_species) :

    delta_z = z_t/np.float(int(n_layers/np.float(fact)))
    h = np.float(fact)*z_t
    n_tropo = int(n_layers/np.float(fact))
    g_t = G*Mp/(R_e+z_t)**2

    alt_array = np.linspace(0,h,n_layers+1)
    alt = np.zeros(n_layers+2)
    alt_mid = np.zeros(n_layers)
    for i_n in range(n_layers) :
        alt_mid[i_n] = (alt_array[i_n+1] +alt_array[i_n])/2.
    alt[0],alt[alt.size-1] = alt_array[0],alt_array[alt_array.size-1]
    alt[1:alt.size-1] = alt_mid

    data = np.zeros((number,n_layers+2))
    g_array = Mp*G/(R_e+alt)**2

    for i_n in range(n_tropo+1) :
        if i_n == 0 :
            r = z_t
        else :
            r = z_t - (i_n-0.5)*delta_z
        data[1,i_n] = T_t +((gamma-1)/gamma)*M*g_t/(R_gp)*(r/(1 - r/(R_e+z_t)))
        #print 'tropo',data[1,i_n]
        data[0,i_n] = P_t*(T_t/data[1,i_n])**(gamma/(1-gamma))
        for i_sp in range(n_species.size) :
            data[2+i_sp,i_n] = x_ratio_species[i_sp]
        data[number-1,i_n] = M

    wh, = np.where(tropo[0] < P_t)
    tropo = tropo[:,wh[:wh.size-1]-1]
    z_tropo = np.zeros(tropo[0].size)
    z_tropo[0] = z_t

    for i_t in range(1,z_tropo.size) :
        if z_tropo[i_t-1] != 0 :
            beta = -R_gp*(tropo[1,i_t]+tropo[1,i_t-1])*(R_e+z_tropo[i_t-1])**2/(2*M*G*Mp)*np.log(tropo[0,i_t]/tropo[0,i_t-1])
            delta = beta/(1-beta/(R_e+z_tropo[i_t-1]))
        if delta > 0 :
            z_tropo[i_t] = z_tropo[i_t-1] + delta

    wh_z, = np.where(z_tropo != 0)
    z_tropo = z_tropo[wh_z]
    tropo = tropo[:,wh_z]

    wh_stra, = np.where(alt_mid > z_tropo[z_tropo.size-1])
    if wh_stra.size != 0 :
        n_strato = wh_stra[0]
        for i_n in range(n_tropo+1,n_strato+1) :
            r = (i_n-0.5)*delta_z
            p_n,i_g,i_c = interpolation(r,z_tropo,tropo[0])
            t_n,i_g,i_c = interpolation(r,z_tropo,tropo[1])
            data[0,i_n] = p_n
            data[1,i_n] = t_n
            #print 'strato', data[1,i_n]
            for i_sp in range(n_species.size) :
                data[2+i_sp,i_n] = x_ratio_species[i_sp]
            data[number-1,i_n] = M

        delta_r = 0.
        z_s = (n_strato-0.5)*delta_z
        for i_n in range(n_strato+1,n_layers+2) :
            data[1,i_n] = data[1,i_n-1]
            #print 'Meso', data[1,i_n]
            if i_n != n_layers+1 :
                delta_r += delta_z
            else :
                delta_r +=delta_z/2.
            data[0,i_n] = data[0,n_strato]*np.exp(-M*g_array[n_strato]/(R_gp*data[1,n_strato])*(delta_r)/(1+delta_r/(R_e+z_s)))
            for i_sp in range(n_species.size) :
                data[2+i_sp,i_n] = x_ratio_species[i_sp]
            data[number-1,i_n] = M
    else :
        for i_n in range(n_tropo+1,n_layers+2) :
            if i_n != n_layers+1 :
                r = (i_n-0.5)*delta_z
            else :
                r = n_layers*delta_z
            p_n,i_g,i_c = interpolation(r,z_tropo,tropo[0])
            t_n,i_g,i_c = interpolation(r,z_tropo,tropo[1])
            data[0,i_n] = p_n
            data[1,i_n] = t_n
            #print 'strato', data[1,i_n]
            for i_sp in range(n_species.size) :
                data[2+i_sp,i_n] = x_ratio_species[i_sp]
            data[number-1,i_n] = M

    H_array = R_gp*data[1]/(M*g_array)
    r_step, x_step, theta_number, reso_alt = delta_z, delta_z, 1, int(h/1000.)
    z_lim = int(alt_array[n_layers]/delta_z)
    z_reso = int(h/delta_z) + 1
    z_array = alt_array

    return data, alt, alt_mid, H_array, g_array, h, delta_z, r_step, x_step, theta_number, reso_alt, z_lim, z_reso, z_array

########################################################################################################################

def tropo_interpolation(R_e,Mp,P_surf,P_h,atmos,M,n_layers,number,n_species,x_ratio_species,extra) :

    press_exp = np.linspace(np.log(P_h), np.log(P_surf), n_layers+1)
    pressure_levels = np.exp(press_exp)[::-1]
    data = np.zeros((number,n_layers+2))
    data[0,0] = P_surf
    data[0,1:n_layers+1] = np.power(10,np.log10(pressure_levels[:-1])+\
                    np.diff(np.log10(pressure_levels))/2.)
    data[0,n_layers+1] = P_h

    if extra.size == 0 :

        wh_atm, = np.where((data[0] >= atmos[0,0,0])*(data[0] <= atmos[0,0,atmos[0,0].size-1]))

        for i_n in wh_atm :
            T_0, i_at, c_at = interpolation(np.log(data[0,i_n]),np.log(atmos[0,0]),atmos[0,1])
            T_1, i_at, c_at = interpolation(np.log(data[0,i_n]),np.log(atmos[1,0]),atmos[1,1])
            data[1,i_n], i_at, c_at = interpolation(np.log(P_surf),np.log(np.array([atmos[0,0,atmos[0,0].size-1],atmos[1,0,atmos[1,0].size-1]])),\
                                    np.array([T_0,T_1]))
        if atmos[1,0,0] < atmos[0,0,0] :
            wh_atm_1, = np.where((data[0] >= atmos[1,0,0])*(data[0] < atmos[0,0,0]))
            for i_n in wh_atm_1 :
                T_1, i_at, c_at = interpolation(np.log(data[0,i_n]),np.log(atmos[1,0]),atmos[1,1])
                data[1,i_n] = T_1
            data[1,wh_atm_1[wh_atm_1.size-1]+1:] = np.ones(n_layers+2-wh_atm_1[wh_atm_1.size-1]-1)*data[1,wh_atm_1[wh_atm_1.size-1]]
        else :
            wh_atm_0, = np.where((data[0] >= atmos[0,0,0])*(data[0] <= atmos[1,0,0]))
            for i_n in wh_atm_0 :
                T_0, i_at, c_at = interpolation(np.log(data[0,i_n]),np.log(atmos[0,0]),atmos[0,1])
                data[1,i_n] = T_0
            data[1,wh_atm_0[wh_atm_0.size-1]+1:] = np.ones(n_layers+2-wh_atm_0[wh_atm_0.size-1]-1)*data[1,wh_atm_0[wh_atm_0.size-1]]

        Ps0, Ps1 = atmos[0,0,atmos[0,0].size-1], atmos[1,0,atmos[1,0].size-1]

        for i_n in range(0,wh_atm[0]) :
            T1, i_t1, c_t1 = interpolation(np.log(data[0,i_n]),np.log(atmos[1,0]),atmos[1,1])
            wh_pente, = np.where(atmos[1,0] > data[0,i_n])
            ok = 0
            i_pente = 0
            while ok == 0 :
                delta_pente = (np.log(atmos[1,0,wh_pente[i_pente]])-np.log(Ps0))/(np.log(Ps1)-np.log(Ps0))
                Ps = np.exp((np.log(data[0,i_n])+(delta_pente-1)*np.log(Ps0))/(delta_pente))
                if Ps <= P_surf or Ps - P_surf < 1.e-9 :
                    ok = 1
                else :
                    i_pente += 1
                Ts, i_ts, c_ts = interpolation(np.log(data[0,i_n]),np.log(np.array([Ps0,atmos[1,0,wh_pente[i_pente]]])),\
                                                np.array([atmos[0,1,atmos[0,0].size-1],atmos[1,1,wh_pente[i_pente]]]))
            data[1,i_n], i_at, c_at = interpolation(np.log(P_surf),np.log(np.array([Ps,Ps1])),np.array([Ts,T1]))

    else :

        P_tropo, T_tropo, gamma = extra[0], extra[1], extra[2]
        wh_atm, = np.where((data[0] >= atmos[0,0])*(data[0] <= P_tropo))

        for i_n in wh_atm :
            data[1,i_n], i_at, c_at = interpolation(np.log(data[0,i_n]),np.log(atmos[0]),atmos[1])

        for i_n in range(0,wh_atm[0]) :
            data[1,i_n] = (P_tropo/data[0,i_n])**((1-gamma)/(gamma))*T_tropo

        data[1,wh_atm[wh_atm.size-1]+1:] = np.ones(n_layers+2-wh_atm[wh_atm.size-1]-1)*data[1,wh_atm[wh_atm.size-1]]

    for i_n in range(n_species.size):
        data[2+i_n] = np.ones(n_layers+2)*x_ratio_species[i_n]
    data[number-1] = np.ones(n_layers+2)*M

    alt_array = np.zeros(n_layers+1)
    alt_array[0] = 0
    g_array = np.zeros(n_layers+1)
    g_array[0] = G*Mp/(R_e**2)
    H_array = np.zeros(n_layers+1)

    H_array[0] = k_B*data[1,0]/(data[number-1,0]*AMU*g_array[0])*1.e-3
    print 'pressure_scaling_taurex says : Scale height : ',R_gp*data[1,0]/(data[number-1,0]*g_array[0]),' m'

    for i_l in range(1,n_layers+1) :
        delta_r = (-1.)*H_array[i_l-1]*np.log(pressure_levels[i_l]/pressure_levels[i_l-1])
        alt_array[i_l] = alt_array[i_l-1] + delta_r
        g_array[i_l] = g_array[0]
        H_array[i_l] = k_B*data[1,i_l]/(data[number-1,i_l]*AMU*g_array[i_l])*1.e-3

    alt_array[n_layers] = alt_array[n_layers-1] + delta_r
    g_array[n_layers] = g_array[0]
    H_array[n_layers] = k_B*data[1,n_layers]/(data[number-1,n_layers]*AMU*g_array[n_layers])*1.e-3

    lim_alt = alt_array[n_layers]
    rupt_alt, h = 0.e+0, lim_alt
    delta_z = h/np.float(n_layers)
    r_step, x_step, theta_number, reso_alt = delta_z, delta_z, 1, int(h/1000)
    z_lim = int(lim_alt/delta_z)
    z_reso = int(h/delta_z) + 1
    z_array = np.arange(h/np.float(delta_z)+1)*float(delta_z)
    alt_mid = np.zeros(n_layers)
    for i_w in range(n_layers) :
        if i_w != n_layers-1 :
            alt_mid[i_w] = (alt_array[i_w+1] + alt_array[i_w])/2.
        else :
            alt_mid[i_w] = alt_array[i_w] + (alt_array[i_w] - alt_array[i_w-1])/2.
    alt = np.zeros(n_layers+2)
    alt[1:n_layers+1] = alt_mid
    alt[0], alt[n_layers+1] = alt_array[0], alt_array[n_layers]

    return data, alt, alt_mid, H_array, g_array, h, delta_z, r_step, x_step, theta_number, reso_alt, z_lim, z_reso, z_array

########################################################################################################################

def tropo_interpolation_generic(R_e,Mp,P_surf,P_h,atmos,n_layers,number,n_species,M_species,c_species,extra) :

    press_exp = np.linspace(np.log(P_h), np.log(P_surf), n_layers+1)
    pressure_levels = np.exp(press_exp)[::-1]
    data = np.zeros((number,n_layers+2))
    data[0,0] = P_surf
    data[0,1:n_layers+1] = np.power(10,np.log10(pressure_levels[:-1])+\
                    np.diff(np.log10(pressure_levels))/2.)
    data[0,n_layers+1] = P_h

    if extra.size == 0 :

        wh_atm, = np.where((data[0] >= atmos[0,0,0])*(data[0] <= atmos[0,0,atmos[0,0].size-1]))

        for i_n in wh_atm :
            T_0, i_at, c_at = interpolation(np.log(data[0,i_n]),np.log(atmos[0,0]),atmos[0,1])
            T_1, i_at, c_at = interpolation(np.log(data[0,i_n]),np.log(atmos[1,0]),atmos[1,1])
            data[1,i_n], i_at, c_at = interpolation(np.log(P_surf),np.log(np.array([atmos[0,0,atmos[0,0].size-1],atmos[1,0,atmos[1,0].size-1]])),\
                                    np.array([T_0,T_1]))
            for i_c in range(c_species.size) :
                c_0, i_at, c_at = interpolation(np.log(data[0,i_n]),np.log(atmos[0,0]),atmos[0,2+i_c])
                c_1, i_at, c_at = interpolation(np.log(data[0,i_n]),np.log(atmos[1,0]),atmos[1,2+i_c])
                data[2+i_c,i_n], i_at, c_at = interpolation(np.log(P_surf),np.log(np.array([atmos[0,0,atmos[0,0].size-1],atmos[1,0,atmos[1,0].size-1]])),\
                                    np.array([c_0,c_1]))
            for i_spe in range(2,n_species.size) :
                spe_0, i_at, c_at = interpolation(np.log(data[0,i_n]),np.log(atmos[0,0]),atmos[0,c_species.size+i_spe])
                spe_1, i_at, c_at = interpolation(np.log(data[0,i_n]),np.log(atmos[1,0]),atmos[1,c_species.size+i_spe])
                data[2+c_species.size+i_spe,i_n], i_at, c_at = interpolation(np.log(P_surf),np.log(np.array([atmos[0,0,atmos[0,0].size-1],atmos[1,0,atmos[1,0].size-1]])),\
                                    np.array([spe_0,spe_1]))
        if atmos[1,0,0] < atmos[0,0,0] :
            wh_atm_1, = np.where((data[0] >= atmos[1,0,0])*(data[0] < atmos[0,0,0]))
            for i_n in wh_atm_1 :
                T_1, i_at, c_at = interpolation(np.log(data[0,i_n]),np.log(atmos[1,0]),atmos[1,1])
                data[1,i_n] = T_1
                for i_c in range(c_species.size) :
                    c_1, i_at, c_at = interpolation(np.log(data[0,i_n]),np.log(atmos[1,0]),atmos[1,2+i_c])
                    data[2+i_c,i_n] = c_1
                for i_spe in range(2,n_species.size) :
                    spe_1, i_at, c_at = interpolation(np.log(data[0,i_n]),np.log(atmos[1,0]),atmos[1,c_species.size+i_spe])
                    data[2+c_species.size+i_spe,i_n] = spe_1
            data[1,wh_atm_1[wh_atm_1.size-1]+1:] = np.ones(n_layers+2-wh_atm_1[wh_atm_1.size-1]-1)*data[1,wh_atm_1[wh_atm_1.size-1]]
            for i_c in range(c_species.size) :
                data[2+i_c,wh_atm_1[wh_atm_1.size-1]+1:] = np.ones(n_layers+2-wh_atm_1[wh_atm_1.size-1]-1)*data[2+i_c,wh_atm_1[wh_atm_1.size-1]]
            for i_spe in range(2,n_species.size) :
                data[2+c_species.size+i_spe,wh_atm_1[wh_atm_1.size-1]+1:] = np.ones(n_layers+2-wh_atm_1[wh_atm_1.size-1]-1)*data[2+c_species.size+i_spe,wh_atm_1[wh_atm_1.size-1]]
        else :
            wh_atm_0, = np.where((data[0] >= atmos[0,0,0])*(data[0] <= atmos[1,0,0]))
            for i_n in wh_atm_0 :
                T_0, i_at, c_at = interpolation(np.log(data[0,i_n]),np.log(atmos[0,0]),atmos[0,1])
                data[1,i_n] = T_0
                for i_c in range(c_species.size) :
                    c_0, i_at, c_at = interpolation(np.log(data[0,i_n]),np.log(atmos[0,0]),atmos[0,2+i_c])
                    data[2+i_c,i_n] = c_0
                for i_spe in range(2,n_species.size) :
                    spe_0, i_at, c_at = interpolation(np.log(data[0,i_n]),np.log(atmos[0,0]),atmos[0,c_species.size+i_spe])
                    data[2+c_species.size+i_spe,i_n] = spe_0
            data[1,wh_atm_0[wh_atm_0.size-1]+1:] = np.ones(n_layers+2-wh_atm_0[wh_atm_0.size-1]-1)*data[1,wh_atm_0[wh_atm_0.size-1]]
            for i_c in range(c_species.size) :
                data[2+i_c,wh_atm_0[wh_atm_0.size-1]+1:] = np.ones(n_layers+2-wh_atm_0[wh_atm_0.size-1]-1)*data[2+i_c,wh_atm_0[wh_atm_0.size-1]]
            for i_spe in range(2,n_species.size) :
                data[2+c_species.size+i_spe,wh_atm_0[wh_atm_0.size-1]+1:] = np.ones(n_layers+2-wh_atm_0[wh_atm_0.size-1]-1)*data[2+c_species.size+i_spe,wh_atm_0[wh_atm_0.size-1]]

        Ps0, Ps1 = atmos[0,0,atmos[0,0].size-1], atmos[1,0,atmos[1,0].size-1]

        for i_n in range(0,wh_atm[0]) :
            T1, i_t1, c_t1 = interpolation(np.log(data[0,i_n]),np.log(atmos[1,0]),atmos[1,1])
            wh_pente, = np.where(atmos[1,0] > data[0,i_n])
            ok = 0
            i_pente = 0
            while ok == 0 :
                delta_pente = (np.log(atmos[1,0,wh_pente[i_pente]])-np.log(Ps0))/(np.log(Ps1)-np.log(Ps0))
                Ps = np.exp((np.log(data[0,i_n])+(delta_pente-1)*np.log(Ps0))/(delta_pente))
                if Ps <= P_surf or Ps - P_surf < 1.e-9 :
                    ok = 1
                else :
                    i_pente += 1
            Ts, i_ts, c_ts = interpolation(np.log(data[0,i_n]),np.log(np.array([Ps0,atmos[1,0,wh_pente[i_pente]]])),\
                                np.array([atmos[0,1,atmos[0,0].size-1],atmos[1,1,wh_pente[i_pente]]]))
            data[1,i_n], i_at, c_at = interpolation(np.log(P_surf),np.log(np.array([Ps,Ps1])),np.array([Ts,T1]))
            for i_c in range(c_species.size) :
                c1, i_t1, c_t1 = interpolation(np.log(data[0,i_n]),np.log(atmos[1,0]),atmos[1,2+i_c])
                cs, i_ts, c_ts = interpolation(np.log(data[0,i_n]),np.log(np.array([Ps0,atmos[1,0,wh_pente[i_pente]]])),\
                                                np.array([atmos[0,2+i_c,atmos[0,0].size-1],atmos[1,2+i_c,wh_pente[i_pente]]]))
                data[2+i_c,i_n], i_at, c_at = interpolation(np.log(P_surf),np.log(np.array([Ps,Ps1])),np.array([cs,c1]))
            for i_spe in range(2,n_species.size) :
                spe1, i_t1, c_t1 = interpolation(np.log(data[0,i_n]),np.log(atmos[1,0]),atmos[1,c_species.size+i_spe])
                spes, i_ts, c_ts = interpolation(np.log(data[0,i_n]),np.log(np.array([Ps0,atmos[1,0,wh_pente[i_pente]]])),\
                                                np.array([atmos[0,c_species.size+i_spe,atmos[0,0].size-1],atmos[1,c_species.size+i_spe,wh_pente[i_pente]]]))
                data[2+c_species.size+i_spe,i_n], i_at, c_at = interpolation(np.log(P_surf),np.log(np.array([Ps,Ps1])),np.array([spes,spe1]))

    else :

        P_tropo, T_tropo, gamma = extra[0], extra[1], extra[2]
        wh_atm, = np.where((data[0] >= atmos[0,0])*(data[0] <= P_tropo))

        for i_n in wh_atm :
            data[1,i_n], i_at, c_at = interpolation(np.log(data[0,i_n]),np.log(atmos[0]),atmos[1])
            for i_c in range(c_species.size) :
                data[2+i_c,i_n], i_at, c_at = interpolation(np.log(data[0,i_n]),np.log(atmos[0]),atmos[2+i_c])
            for i_spe in range(2,n_species.size) :
                data[2+c_species.size+i_spe,i_n], i_at, c_at = interpolation(np.log(data[0,i_n]),np.log(atmos[0]),atmos[c_species.size+i_spe])

        for i_n in range(0,wh_atm[0]) :
            data[1,i_n] = (P_tropo/data[0,i_n])**((1-gamma)/(gamma))*T_tropo
        for i_c in range(c_species.size) :
            data[2+i_c,0:wh_atm[0]] = np.ones(wh_atm[0])*data[2+i_c,wh_atm[0]]
        for i_spe in range(2,n_species.size) :
            data[2+c_species.size+i_spe,0:wh_atm[0]] = np.ones(wh_atm[0])*data[c_species.size+i_spe,wh_atm[0]]

        data[1:,wh_atm[wh_atm.size-1]+1:] = np.ones((number-1,n_layers+2-wh_atm[wh_atm.size-1]-1))*data[1:,wh_atm[wh_atm.size-1]]

    for i_n in range(n_layers+2):
        data[2+c_species.size,i_n] = (1.-np.nansum(data[4+c_species.size:4+c_species.size+n_species.size,i_n]))/(1.+ratio_HeH2)
        data[3+c_species.size,i_n] = data[2+c_species.size,i_n]*ratio_HeH2
        data[number-1] = np.nansum(data[2+c_species.size:2+c_species.size+n_species.size,i_n]*M_species)

    alt_array = np.zeros(n_layers+1)
    alt_array[0] = 0
    g_array = np.zeros(n_layers+1)
    g_array[0] = G*Mp/(R_e**2)
    H_array = np.zeros(n_layers+1)

    H_array[0] = k_B*data[1,0]/(data[number-1,0]*AMU*g_array[0])*1.e-3
    print 'pressure_scaling_taurex says : Scale height : ',R_gp*data[1,0]/(data[number-1,0]*g_array[0]),' m'

    for i_l in range(1,n_layers+1) :
        delta_r = (-1.)*H_array[i_l-1]*np.log(pressure_levels[i_l]/pressure_levels[i_l-1])
        alt_array[i_l] = alt_array[i_l-1] + delta_r
        g_array[i_l] = g_array[0]
        H_array[i_l] = k_B*data[1,i_l]/(data[number-1,i_l]*AMU*g_array[i_l])*1.e-3

    alt_array[n_layers] = alt_array[n_layers-1] + delta_r
    g_array[n_layers] = g_array[0]
    H_array[n_layers] = k_B*data[1,n_layers]/(data[number-1,n_layers]*AMU*g_array[n_layers])*1.e-3

    lim_alt = alt_array[n_layers]
    rupt_alt, h = 0.e+0, lim_alt
    delta_z = h/np.float(n_layers)
    r_step, x_step, theta_number, reso_alt = delta_z, delta_z, 1, int(h/1000)
    z_lim = int(lim_alt/delta_z)
    z_reso = int(h/delta_z) + 1
    z_array = np.arange(h/np.float(delta_z)+1)*float(delta_z)
    alt_mid = np.zeros(n_layers)
    for i_w in range(n_layers) :
        if i_w != n_layers-1 :
            alt_mid[i_w] = (alt_array[i_w+1] + alt_array[i_w])/2.
        else :
            alt_mid[i_w] = alt_array[i_w] + (alt_array[i_w] - alt_array[i_w-1])/2.
    alt = np.zeros(n_layers+2)
    alt[1:n_layers+1] = alt_mid
    alt[0], alt[n_layers+1] = alt_array[0], alt_array[n_layers]

    return data, alt, alt_mid, H_array, g_array, h, delta_z, r_step, x_step, theta_number, reso_alt, z_lim, z_reso, z_array