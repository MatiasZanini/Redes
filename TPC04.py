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
import matplotlib.cm as cm
import pandas as pd
import igraph as ig
import networkx as nx

#%%

'''
Punto 3)
'''

'''
Item a)
'''

# -------------- Carga de datos-------------------

path = 'C:/Users/Mati/Documents/GitHub/Redes/TC04_Data/'

file = 'geneX.csv'

genes_df = pd.read_csv(path+file) # Abrimos el archivo como objeto de pandas

genes_df = genes_df.T # Trasponemos la matriz para que las columnas contenga la evolucion temporal de cada proteina

protein_names = genes_df.iloc[0]

genes_df = pd.DataFrame(genes_df.values[1:], columns = protein_names) # Creamos el dataframe de pandas con los datos

genes_df= genes_df.astype(float)  # Convertimos los elementos del dataframe en punto flotante

'''
NOTA: las filas no se encuentran etiquetadas. Los índices representan el paso temporal de 0 a 11. Las columnas tienen
la etiqueta de la proteína en cuestión.
'''

#%%

correlacion = genes_df.corr(method='pearson', min_periods=1) # Calcula la correlacion de Pearson para las columnas

correlacion.set_index([correlacion.iloc[0], correlacion.columns[0]]) # Agregamos nombres a las columnas y filas

'''
NOTA: la matriz de correlación tiene unos en su diagonal. Esto tiene sentido ya que la correlación de Pearson entre dos
elementos iguales (elementos en la diagonal, por ejemplo), es 1.
'''

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

red_coexp = nx.convert_matrix.from_pandas_adjacency(coexp) # Convertimos la matriz de adyacencia en un grafo de networkx

nx.draw_kamada_kawai(red_coexp) # Graficamos la red para visualizarla

coexp_numpy = coexp.values

red_coexp_ig = ig.Graph.TupleList(red_coexp.edges(), directed=False) # Creamos la misma red, pero como objeto de Igraph

# layout = red_coexp_ig.layout_kamada_kawai()
# layout = red_coexp_ig.layout("kamada_kawai")
# ig.plot(red_coexp_ig, layout= layout)


#%%

'''
Item c)
'''
# -----------Funciones utiles---------------------


def ig_part2dict(Red_igraph, particion_igraph):
    
    '''
    Convierte una partición en comunidades a un diccionario
    '''
    
    particion_dict = {}
    
    for cluster in range(len(particion_igraph)):
    
        for nodo in Red_igraph.vs(particion_igraph[cluster])['name']:
        
            particion_dict.update({nodo:cluster})
  
    return particion_dict



def graficar_particion(Red, particion_diccionario, posiciones, tamano_nodo = 200, colormap = 'viridis'):
    
    '''
    Grafica una red de networkx con colores segun su partició dada por un diccionario según un dado set de posiciones 
    para los nodos.
    '''
    
    cmap = cm.get_cmap(colormap, max(particion_diccionario.values())+1) # viridis es el mapa de colores
    
    # grafico los nodos
    nx.draw_networkx_nodes(Red, posiciones, particion_diccionario.keys(), node_size = tamano_nodo,
                           cmap=cmap, node_color = list(particion_diccionario.values()), with_labels = False)
    
    # grafico los enlaces aparte
    nx.draw_networkx_edges(Red, pos = posiciones, alpha=0.5)


#%%

# Infomap:

infomap_ig = red_coexp_ig.community_infomap() # Creamos la partición mediante el algoritmo de Infomap

infomap_dict = ig_part2dict(red_coexp_ig, infomap_ig) # Convertimos la partición por Infomap a un diccionario


# Fast Greedy:

fg_dendograma = red_coexp_ig.community_fastgreedy(weights=None) # Dendograma con el FastGreedy que maximiza Q

fg_ig = fg_dendograma.as_clustering() # Comunidades dadas por Fast Greedy como objeto de Igraph

fg_dict = ig_part2dict(red_coexp_ig, fg_ig)


'''
Estimación de las modularidades:
'''

Q_infomap = red_coexp_ig.modularity(infomap_ig)

Q_fg = red_coexp_ig.modularity(fg_ig)

print('La modularidad por infomap es:', Q_infomap)

print('La modularidad por Fast Greedy es:', Q_fg)

#%%

#---------------- Graficamos las particiones --------------------

posiciones = nx.kamada_kawai_layout(red_coexp)

colormap = 'Accent'

tamano_nodo = 35

plt.figure(figsize=(15,15))

plt.subplot(1,2,1)

plt.title('Coexpresión Génica por Infomap')

graficar_particion(red_coexp, infomap_dict, posiciones, tamano_nodo, colormap)

plt.subplot(1,2,2)

plt.title('Coexpresión Génica por Fast Greedy')

graficar_particion(red_coexp, fg_dict, posiciones, tamano_nodo, colormap)

'''
Viendo los grafos graficados, puede observarse que el mecanismo de Fast Greedy pareciera tener una menor granularidad
en la distribución de sus comunidades. Esto quiere der, que se observa una menor cantidad de particiones y de mayor
tamaño. Para visualizar mejor esto, realizaremos un histograma mostrando la cantidad de proteinas incluida en cada
partición para cada uno de los dos métodos.
'''

#%%

conteo_fg = list(fg_dict.values())

cant_comunidades_fg = len(fg_ig)

bins_fg = np.linspace(0, cant_comunidades_fg-1, cant_comunidades_fg)

conteo_infomap = list(infomap_dict.values())

cant_comunidades_infomap = len(infomap_ig)

bins_infomap = np.linspace(0, cant_comunidades_infomap-1, cant_comunidades_infomap)

plt.figure()

plt.hist(conteo_fg, bins= bins_fg, density=False, facecolor='blue', alpha=0.9, ec='black', label = 'Fast Greedy')

plt.hist(conteo_infomap, bins= bins_infomap, density=False, facecolor='red', alpha=0.4, ec='black', label = 'Infomap')

plt.xlabel('Comunidades')

plt.ylabel('Cantidad de nodos')

plt.legend()

plt.grid()

print('La partición por Fast Greedy contiene', cant_comunidades_fg, 'comunidades.')

print('La partición por Infomap contiene', cant_comunidades_infomap, 'comunidades.')

'''            
Efectivamente, tal como se observa en el histograma, la partición realizada por Fast Greedy, posee una granularidad 
menor. Esto se traduce en una menor cantidad total de comunidades, y una mayor población en las mismas.
'''
