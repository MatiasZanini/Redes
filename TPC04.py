# -*- coding: utf-8 -*-
"""
                **Introducción a Redes Complejas en Biología de Sistemas**
                        Trabajo Computacional 4 (Entrega 17/06)


Grupo: Camila Sanz y Matías Zanini.
"""
################################################################################
#                                 PAQUETES 
################################################################################
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#%%

'''
Punto 3)
'''


# Item a)


#%% -------------- Carga de datos

path = 'C:/Users/Mati/Documents/GitHub/Redes/TC04_Data/'

file = 'geneX.csv'

genes_df = pd.read_csv(path+file) # Abrimos el archivo como objeto de pandas

genes_df = genes_df.T # Trasponemos la matriz para que las columnas contenga la evolucion temporal de cada proteina

protein_names = genes_df.iloc[0]

genes_df = pd.DataFrame(genes_df.values[1:], columns = protein_names) # Creamos el dataframe de pandas con los datos

# NOTA: las filas no se encuentran etiquetadas. Los índices representan el paso temporal de 0 a 11. Las columnas tienen la
# etiqueta de la proteina en cuestion.


#%%