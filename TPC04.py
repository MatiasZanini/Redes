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

'''
Item a)
'''

#%% -------------- Carga de datos

path = 'C:/Users/Mati/Documents/GitHub/Redes/TC04_Data/'

file = 'geneX.csv'

genes_df = pd.read_csv(path+file) # Abrimos el archivo como objeto de pandas

genes_df = genes_df.T # Trasponemos la matriz para que las columnas contenga la evolucion temporal de cada proteina

protein_names = genes_df.iloc[0]

genes_df = pd.DataFrame(genes_df.values[1:], columns = protein_names) # Creamos el dataframe de pandas con los datos

genes_df= genes_df.astype(float)  # Convertimos los elementos del dataframe en punto flotante

# NOTA: las filas no se encuentran etiquetadas. Los índices representan el paso temporal de 0 a 11. Las columnas tienen la
# etiqueta de la proteína en cuestión.


#%%

correlacion = genes_df.corr(method='pearson', min_periods=1) # Calcula la correlacion de Pearson para las columnas

correlacion.set_index([correlacion.iloc[0], correlacion.columns[0]]) # Agregamos nombres a las columnas y filas

# NOTA: la matriz de correlación tiene unos en su diagonal. Esto tiene sentido ya que la correlación de Pearson entre dos
# elementos iguales (elementos en la diagonal, por ejemplo), es 1.

similaridad = (1 + correlacion)/2 # Matriz de similariadad

'''
La matriz de similaridad se define de esta manera ya que, según la correlación de Pearson, la máxima correlación
posible entre dos elementos se representa con un 1. El valor 0 indica que la correlación es inexistente y un valor -1
indica una correlacion inversa entre los elementos. La matriz de similitud normaliza esta noción de correlacion.
De esta manera, el máximo valor de correlación se alcanzará en 1 (elementos completamente iguales) mientras que la 
correlación inversa, indica un grado mínimo de similaridad, reperesentado por el valor nulo.
'''

#%%


'''
item b)
'''

coexp = pd.DataFrame(0, index = protein_names, columns= protein_names) # Inicializamos un dataframe lleno de ceros

for i in protein_names:
    
    for j in protein_names:
        
        if similaridad[i][j] >= 0.95: # Creamos la matriz de coexperesion genica
            
            coexp[i][j] = 1


'''
Queda: 
    
    1. crear la particion por fastgreedy (facil con igraph en el colab) 
    
    2. idem por infomap
    
    3. calcular la modularidad de cada particion  (sacar del colab)
    
    4. graficar las redes con las comunidades coloereadas (sacar del colab)
    
    5. preguntar que es granularidad y tipo para una comunidad
    
'''


'''            
Calcule la partición en clusters de dicha red mediante los métodos infomap y fastgreedy. Estime
la modularidad de ambas particiones. Visualice ambas redes y compare el tipo y granularidad de
las particiones obtenidas.
'''




# dendograma_fast_greedy = Red_delfines_igraph.community_fastgreedy(weights=None)
# print(dendograma_fast_greedy)
# print('La partición de Modularidad optima tiene '+str(dendograma_fast_greedy.optimal_count)+' comunidades')
# print('Así vemos las comunidades:')
# print(dendograma_fast_greedy.as_clustering())


















