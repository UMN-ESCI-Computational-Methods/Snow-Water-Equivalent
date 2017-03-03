
#Keenen Francois-King
#Assignment 1

#import modules
import numpy as np
from scipy.optimize import fsolve

#constants
rho_w = 1000
rho_ds = 370
rho_i = 917
freq = 1.57542e9
eps_a_prime = 1.0
eps_i_prime = 3.18
eps_w_prime = 88
eps_w_dub_prime = 9.8
eps_0 = 8.8541878176e-12
mu_0 = 1.25663706e-6
Zv = np.sqrt(mu_0/eps_0)

#user input from data that is collected at the site of measurement
Im1 = float(input("What is the normalized intensity measured at GPS1 above the snowpack (dB)?\n"))
Im23 = float(input("What is the normalized mean intensity of GPS2 and GPS3 below the snowpack (dB)?\n"))
d = float(input("What is the snow depth? (meters)\n"))
theta_elev = float(input("What is the elevation angle received at GPS1 from the satellite? (degrees)\n"))
print("\n")

#converting theta_elev to theta_0
theta_0 = float(np.radians(90-theta_elev))


#Function for calculating LWE using the empirical formula for eps_s_prime via Sihvola and Tiuri
def f_sihvola(theta_w_sihvola):
    
    eps_s_prime = 1 + (1.7e-3)*rho_ds+(7.0e-7)*rho_ds**2+(8.7e-2)*theta_w_sihvola+(7.0e-3)*theta_w_sihvola**2
    eps_s_dub_prime = (freq/(10**9))*((1.0e-3)*theta_w_sihvola+((8.0e-5)*theta_w_sihvola**2))*eps_w_dub_prime
    ds = (d)/(np.cos(np.arcsin(np.sin(theta_0/np.sqrt(eps_s_prime)))))
    theta_refr = (np.arcsin(np.sin(theta_0/np.sqrt(eps_s_prime))))
    Zs = np.sqrt(mu_0/(eps_0*(eps_s_prime+eps_s_dub_prime)))
    r_parr = (Zv*np.cos(theta_0)-Zs*np.cos(theta_refr))/(Zv*np.cos(theta_0)+Zs*np.cos(theta_refr))
    r_perp = (Zs*np.cos(theta_0) - Zv*np.cos(theta_refr))/(Zs*np.cos(theta_0) + Zv*np.cos(theta_refr))
    Ir = ((r_perp**2 + r_parr**2)/2)*Im1

    attenuation1 = np.sqrt(mu_0/(eps_s_prime*eps_0))*eps_s_dub_prime*eps_0*2*np.pi*freq
    attenuation2 = -np.log(Im23/(Im1-Ir))/ds
    set_equal = attenuation1 - attenuation2
    
    return set_equal

theta_w_sihvola = fsolve(f_sihvola, 0.01)
theta_w_sihvola, f_sihvola(theta_w_sihvola)


#Function for calculating LWE using the empirical formula for eps_s_prime via Denoth
def f_denoth(theta_w_denoth):
    
    rho_ws = rho_ds +(0.01*theta_w_denoth*rho_w)
    eps_s_prime = 1+(1.92e-3)*rho_ws+(4.4e-7)*(rho_ws)**2 + (1.87e-1)*theta_w_denoth+(4.5e-3)*theta_w_denoth**2
    eps_s_dub_prime = (freq/(10**9))*((1.0e-3)*theta_w_sihvola+((8.0e-5)*theta_w_sihvola**2))*eps_w_dub_prime
    ds = (d)/(np.cos(np.arcsin(np.sin(theta_0/np.sqrt(eps_s_prime)))))
    theta_refr = (np.arcsin(np.sin(theta_0/np.sqrt(eps_s_prime))))
    Zs = np.sqrt(mu_0/eps_0*(eps_s_prime+eps_s_dub_prime))
    r_parr = (Zv*np.cos(theta_0)-Zs*np.cos(theta_refr))/(Zv*np.cos(theta_0)+Zs*np.cos(theta_refr))
    r_perp = (Zs*np.cos(theta_0) - Zv*np.cos(theta_refr))/(Zs*np.cos(theta_0) + Zv*np.cos(theta_refr))
    Ir = ((r_perp**2 + r_parr**2)/2)*Im1

    attenuation1 = np.sqrt(mu_0/(eps_s_prime*eps_0))*eps_s_dub_prime*eps_0*2*np.pi*freq
    attenuation2 = -np.log(Im23/(Im1-Ir))/ds
    set_equal = attenuation1 - attenuation2
    
    return set_equal
    
theta_w_denoth = fsolve(f_denoth, 0.01)
theta_w_denoth, f_denoth(theta_w_denoth)   

 
#Function for calculating LWE using the empirical formula for eps_s_prime via Roth et al.
def f_roth(theta_w_roth):
    
    eps_s_prime = (0.01*theta_w_roth*(eps_w_prime**0.5)+(rho_ds/rho_i)*(eps_i_prime**0.5)+(1-(rho_ds/rho_i)-0.01*theta_w_roth)*(eps_a_prime**0.5))**2
    eps_s_dub_prime = (freq/(10**9))*((1.0e-3)*theta_w_sihvola+((8.0e-5)*theta_w_sihvola**2))*eps_w_dub_prime
    ds = (d)/(np.cos(np.arcsin(np.sin(theta_0/np.sqrt(eps_s_prime)))))
    theta_refr = (np.arcsin(np.sin(theta_0/np.sqrt(eps_s_prime))))
    Zs = np.sqrt(mu_0/eps_0*(eps_s_prime+eps_s_dub_prime))
    r_parr = (Zv*np.cos(theta_0)-Zs*np.cos(theta_refr))/(Zv*np.cos(theta_0)+Zs*np.cos(theta_refr))
    r_perp = (Zs*np.cos(theta_0) - Zv*np.cos(theta_refr))/(Zs*np.cos(theta_0) + Zv*np.cos(theta_refr))
    Ir = ((r_perp**2 + r_parr**2)/2)*Im1

    attenuation1 = np.sqrt(mu_0/(eps_s_prime*eps_0))*eps_s_dub_prime*eps_0*2*np.pi*freq
    attenuation2 = -np.log(Im23/(Im1-Ir))/ds
    set_equal = attenuation1 - attenuation2
    
    return set_equal

theta_w_roth = fsolve(f_roth, 0.01)
theta_w_roth, f_roth(theta_w_roth)

print("Volumetric LWC for Sihvola and Tiuri (percent volume): {}%".format("%.5f" % theta_w_sihvola[0]),'\n')
print("Volumetric LWC for Denoth (percent volume): {}%".format("%.5f" % theta_w_denoth[0]),'\n')
print("Volumetric LWC for Roth (percent volume): {}%".format("%.5f" % theta_w_roth[0]),'\n')