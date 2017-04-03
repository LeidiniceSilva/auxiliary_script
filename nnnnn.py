import matplotlib as mpl ; mpl.use('Agg')  # Descomente para não mostrar a janela em cada plot
import numpy as np
from numpy import loadtxt
import os
import matplotlib.pyplot as plt

# Fazendo Dezembro

xaxis = ['Jan', 'Fev', 'Mar', 'Abr', 'Mai', 'Jun', 'Jul', 'Ago', 'Set', 'Out', 'Nov', 'Dez']

path_esi = "/home/junior/funceme/projetos_aux/esi/data/asc_mon/"
path_pet = "/home/junior/funceme/projetos_aux/esi/pet/"
path_out = "/home/junior/Desktop/"
pet_var = loadtxt(path_pet + 'bac_value_pet_012016_122016.asc', delimiter=";")
pet_bacs = np.genfromtxt(path_pet + 'bac_names_pet.asc',dtype='str')

filesesi = os.listdir(path_esi)
t = np.arange(0,12)

for i,bac in enumerate(pet_bacs):

    try:

        esifile = [ f for f in filesesi if f.find(bac) != -1 ][0]
        esi = loadtxt(path_esi + esifile, delimiter=",")[-12:,2]
        #esi = loadtxt(path_esi + esifile, delimiter=",")[10:22,2] #refer out
        pet = pet_var[i,:]

        fig, ax1 = plt.subplots(figsize=(120,56))
        plt.xticks(t+0.3, xaxis,fontsize=194,fontweight='bold')
        plt.yticks(fontweight='bold', fontsize=194)

        plt.rc('font', weight='bold')
        ln1 =  ax1.plot(t+0.3, pet, 'b-', label="ETP", linewidth=20)
        ax1.set_xlabel(u'Meses', fontsize=200, fontweight='bold')
        ax1.set_ylabel(u'ETP (mm)', fontsize=200, fontweight='bold')
        ax1.set_ylim(0, 450)

        ax2 = ax1.twinx()
        ln2 = ax2.plot(t+0.3, esi, 'r-', label="ESI", linewidth=20)
        ax2.set_ylabel('ESI ($\sigma$)', fontsize=200, fontweight='bold')
        ax2.set_ylim(-3, 3)
        lns = ln1+ln2
        labs = [l.get_label() for l in lns]
        ax1.legend(lns, labs, loc="best", frameon=False, fontsize=200)
        ax1.tick_params(axis='both', which='major', labelsize=194,length=30, width=10, pad=80)
        ax2.tick_params(axis='both', which='major', labelsize=194, pad=30)
        figname = "{0}{1}.png".format(path_out, bac)
        plt.savefig(figname, bbox_inches='tight', dpi=13)
        plt.close()

    except:
        print u"Bacia {0} não encontrada!!!".format(bac)
        exit()
