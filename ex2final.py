
git clone https://gitlab.developers.cam.ac.uk/ch/thom/part2programming.git
import matplotlib.pyplot as plt
import numpy as np
import glob
from scipy.optimize import curve_fit
from scipy.constants import speed_of_light, atomic_mass


line_angle_length = 98
Hartree_to_Joule = 4.35974 * 10**(-18)
Angstrom_to_metre = 10**(-10)
Degree_to_rad = (np.pi)/180
freq_to_wn = 1 / (speed_of_light * 10 ** 2)

def geometry_extractor(line):
    """Returns the specified bond length and angle for a triatomic for this type of output file."""
    data = line.split()
    
    length = float(data[2])
    angle = float(data[4])

    return (length, angle)

def energy_extractor(line):
    """Returns the specified energy for a triatomic for this type of output file."""
    data = line.split()
    energy = float(data[4])

    return (energy)

def extract_data(data):
    """Returns bond lengths, bond angles and energies for each calculation.
    
    Keyword arguments:
    data -- list of output files for a given molecule.

    Return values:
    h2xbondlengths_ret -- list of bond lengths extracted (Angstroms).
    h2xbondangles_ret -- list of bond angles extracted (Degrees).
    h2xenergies_ret -- list of total energies extracted (Hartree).
    index_min -- index of minimum energy value.
    """

    h2xbondlengths = []
    h2xbondangles = []
    h2xenergies = []

    for file in data:
        f = open(file, "r")
        i=0
        geom_energ=[]
        for line in f:
            i+=1
            if i == line_angle_length:
                length, angle = geometry_extractor(line)
                geom_energ.append(length)
                geom_energ.append(angle)
            if "SCF Done" in line:
                energy = energy_extractor(line)
                geom_energ.append(energy)

                

        h2xbondlengths.append(geom_energ[0])
        h2xbondangles.append(geom_energ[1])
        h2xenergies.append(geom_energ[2])

    index_min = np.argmin(h2xenergies)
    
    h2xbondlengths_ret = np.array(h2xbondlengths)
    h2xbondangles_ret = np.array(h2xbondangles)
    h2xenergies_ret = np.array(h2xenergies)

    return h2xbondlengths_ret, h2xbondangles_ret, h2xenergies_ret, index_min


def fit_data_selection(lengths, angles, energies):
    """Returns a selection of data around the minimum.
    
    Keyword arguments:
    lengths -- list of bond lengths extracted (Angstroms).
    angles -- list of bond angles extracted (Degrees).
    energies -- list of total energies extracted (Hartree).

    Return values:
    selected_geom -- list of lengths (m) and angles (rad) in specified range.
    selected_e -- list of energies corresponding to the geometries above (J).
    """
    selected_l=[]
    selected_a=[]
    selected_e=[]


    for i in range(len(energies)):
        
        if np.absolute(angles[i] - angles[index_min]) < 6 and \
            np.absolute(lengths[i] - lengths[index_min]) < 0.06:

            selected_l.append(lengths[i])
            selected_a.append(angles[i])
            selected_e.append(energies[i])

    selected_geom = [np.array(selected_l) * Angstrom_to_metre, np.array(selected_a) * Degree_to_rad]
    selected_e_SI = np.array(selected_e) * Hartree_to_Joule

    return selected_geom, selected_e_SI 
    
def harm_oscill(geom, k_r, k_theta, length_min, angle_min, E_min):
    """Harmonic Oscillator Function to fit against."""
    length = geom[0]
    angle = geom[1]
    return (((k_r * (length - length_min) ** 2) + (k_theta * (angle - angle_min) ** 2 )) / 2 + E_min)

def freq_calculator(k_opt, length_eq):
    """Calculate and print stretching frequencies"""
    k_r = k_opt[0]
    k_theta = k_opt[1]

    v_stretch = 1 / (2 * np.pi) * np.sqrt(k_r / (2 * atomic_mass))
    print(f"The frequency of the symmetric stretch mode is {v_stretch * freq_to_wn} cm-1.")

    v_bend = 1 / (2 * np.pi) * np.sqrt(k_theta / (0.5 * atomic_mass * length_eq ** 2))
    print (f"The frequency of the bending mode is {v_bend * freq_to_wn} cm-1.")


directories = glob.glob("./*outfiles")

ValidInput = False

system_input = input(f"""Please choose a triatomic molecule: 
            1: {directories[0][2:5]}
            2: {directories[1][2:5]}
            """)

if (system_input == "1") or (system_input == "2"):
    ValidInput = True

system_id = int(system_input)

if ValidInput is True:

    chosen_system = directories[system_id - 1]
    output = glob.glob(f"{chosen_system}/*")
    output.sort()

else: print("Invalid Input. Please restart program and try again.")

lengths, angles, energies, index_min = extract_data(output)

print(f"The equilibrium energy is: {energies[index_min]} at a bond length of {lengths[index_min]} Angstrom and a bond angle of {angles[index_min]} degrees")


lengths_SI = lengths * Angstrom_to_metre
angles_SI = angles * Degree_to_rad
energies_SI = energies * Hartree_to_Joule

print(f"The equilibrium energy is: {energies_SI[index_min]} at a bond length of {lengths_SI[index_min]} Angstrom \
      and a bond angle of {angles_SI[index_min]} degrees")


fit_geom, fit_energies = fit_data_selection(lengths, angles, energies)

print(fit_geom, fit_energies)



popt, pcov = curve_fit(lambda geom, k_r, k_theta: 
            harm_oscill(geom ,k_r, k_theta, lengths_SI[index_min], angles_SI[index_min], energies_SI[index_min]),
            fit_geom, fit_energies)
freq_calculator(popt, lengths_SI[index_min])

molecule = directories[system_id - 1][2:5]

fig = plt.figure(figsize=(12, 9))
ax = fig.add_subplot(111, projection='3d')

surf = ax.plot_trisurf(angles, lengths, energies, cmap='coolwarm', linewidth=0.1, antialiased=True, edgecolor='gray')


ax.set_xlabel('Angle / Degrees')
ax.set_ylabel('O-H bond length / Angstroms ')
ax.set_zlabel('Energy / Hartrees')
ax.set_title(f'{molecule} Potential Energy Surface')

ax.view_init(elev=29, azim=-200)

plt.savefig(f"./images/{molecule}.png", bbox_inches = "tight")

