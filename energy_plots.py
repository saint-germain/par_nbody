import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib.ticker as mtick
from tqdm import tqdm
import sys, os

# Usar formato .tex
#rc('text', usetex=True)
#rc('font', family='Computer Modern Roman')

# Accendiendo a la carpeta con los datos
folder_data = sys.argv[1]
systems_folder = os.listdir(folder_data)

folder_save = sys.argv[2]

# Creando carpeta donde se guardan las graficas
try:
    os.mkdir(folder_save)
except:
    print('Carpeta ya creada')

# Parametros
tf = 1e6
timeline = np.linspace(0, tf, 10000)

G = 6.674e-11

# Buscando energias
for folder in tqdm(systems_folder):
    path_txt = os.path.join(folder_data, folder)
    
    txt_nmass = np.genfromtxt(str(path_txt) + '/Nmass.txt', delimiter=',')
    mstar = txt_nmass[0][8]

    mu_paramater = G*mstar

    orbital_elements = [orb_ele for orb_ele in os.listdir(os.path.join(folder_data, folder)) if orb_ele.startswith('OrbitalElements')]

    temp_energy = np.array([])


    for energy_index in orbital_elements:
        txt_orb_elements = np.genfromtxt(path_txt + '/' + energy_index, delimiter=',')
        semimajor = txt_orb_elements[:,2]


        specific_orbital_energy = -mu_paramater/(2*semimajor)
        temp_energy = np.append(temp_energy, np.sum(specific_orbital_energy))


    # Energia inicial
    energy_ini = temp_energy[0]

    # Graficando
    fig = plt.figure(figsize=(12, 5))
    ax = fig.add_subplot(111)

    # ax.set_title('Excentricidad media vs Tiempo')
    # ax.ticklabel_format(axis='both', style='sci', useMathText=True)
    
    ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))

    ax.grid(ls='--', alpha=0.3)
    # ax.set_xscale('log')
    # ax.set_yscale('log')
    
    # ax_ecc_plot.set_ylim(0, Mi.max() + 5)
    
    ax.set_xlabel('Tiempo (Años)')
    ax.set_ylabel('Error energía relativa')
    # ax.plot(timeline[:-1], temp_energy)
    try:
        ax.plot(timeline, (temp_energy - energy_ini)/energy_ini)
    except:
        print('Dimensiones no compatibles')
    fmt= 'svg' # Formato .svg
    fname=folder_save + '/' + '{}_energy.{}'.format(folder, fmt)
    fig.savefig(fname, format=fmt, bbox_inches='tight')
    plt.close(fig)
