#Keenen Francois-King
#Project 1



#user input from data that is collected at the site of measurement
Im1 = float(input("What is the normalized intensity measured at GPS1 above the snowpack (dB)?\n"))
Im23 = float(input("What is the normalized mean intensity of GPS2 and GPS3 below the snowpack (dB)?\n"))
d = float(input("What is the snow depth? (meters)\n"))
theta_elev = float(input("What is the elevation angle received at GPS1 from the satellite? (degrees)\n"))
print("\n")

import liquid_water_content

solve = liquid_water_content.LWE(Im1, Im23, d, theta_elev)
solve.solve_equations()







