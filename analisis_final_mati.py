# -*- coding: utf-8 -*-
"""
Created on Thu Jul  9 16:54:37 2020

@author: Mati
"""

import os
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import networkx as nx
from networkx.algorithms import bipartite
import igraph as ig
import matplotlib.cm as cm

#%%-------Unificar dataframes---------------

# Categorías: [humor negro, humor verde, humor de serie, humor interno, humor de identificacion] (del 1 al 5)

save_path = 'C:/Users/Mati/Documents/GitHub/Redes/datos/data_frames/'

df1 = pd.read_pickle(save_path + 'primera_tanda.p')

df2 = pd.read_pickle(save_path + 'segunda_tanda.p')

df3 = pd.read_pickle(save_path + 'datos_red.p')

borrables = [3,8,9,20,47,46, 48, 49, 50, 52, 54,56,53]

df3.drop(borrables, inplace = True) # Eliminamos posts mal categorizados.

df1.replace({'5' : '3', '6': '4', '7':'5'}, inplace=True)

df2.replace({'5' : '3', '6': '4', '7':'5'}, inplace=True)

data_total = df1

data_total = data_total.append(df2, ignore_index=True)

data_total = data_total.append(df3, ignore_index=True)

data_total.drop_duplicates(subset ='url', keep = 'first', inplace = True)

data_total.to_pickle(save_path + 'data_total.p')

#data_ordenada = data_total.sort_values('categoria')


#%%

save_path = 'C:/Users/Mati/Documents/GitHub/Redes/datos/data_frames/'

save_path_red = 'C:/Users/Mati/Documents/GitHub/Redes/datos/redes_gml/'

save_name = 'data_total.p'


data = pd.read_pickle(save_path+save_name)

enfermitos = []

posts = []

for i in range(len(data)):
    
    posts.append(data['url'][i])
    
    poster = data['poster'][i]
    
    reacters = list(data['reacters'][i])
    
    for j in reacters:
        
        if j not in enfermitos:
            
            enfermitos.append(j)
            
    if poster not in enfermitos:
        
        enfermitos.append(poster)


red_completa = nx.DiGraph()

red_completa.add_nodes_from(enfermitos, bipartite = 0)

red_completa.add_nodes_from(posts, bipartite = 1)

for i in range(len(data)):
    
    post = data['url'][i]
    
    poster = data['poster'][i]
    
    reacters = list(data['reacters'][i])
    
    red_completa.add_edge(poster, post)

    for j in reacters:
        
        red_completa.add_edge(post, j)


nx.write_gml(red_completa, save_path_red+'red_completa.gml')


#%% 

############
#ANALISIS
############

# Carga de datos

save_path_red = 'C:/Users/Mati/Documents/GitHub/Redes/datos/redes_gml/'

red_name = 'red_completa.gml'

red_completa = nx.read_gml(save_path_red + red_name)

#%%-----------Funciones---------------


def nodo_conectado(red, nodo1, nodo2):
    
    return nodo1 in red.neighbors(nodo2) # Devuelve True si estan conectados o False sino


def proyect_sin_posters(red):
    
    proyectada = nx.Graph()
    
    enfermitos, posts = nx.bipartite.sets(red)
    
    for post in posts:
        
        for reacter in red.neighbors(post):
            
            proyectada.add_node(reacter, bipartite=0)
            
            proyectada.add_node(post, bipartite=1)
            
            proyectada.add_edge(reacter, post)
    
    
    
    return proyectada
    

def proyect_enfermitos(red, tol):
    
    '''
    si 2 enfermitos le dieron like a "tol" o mas posts iguales --> enlace
    '''
    
    enfermitos, posts = nx.bipartite.sets(red)
    
    red_enfermita = bipartite.weighted_projected_graph(red, enfermitos)
    
    enlaces = red_enfermita.edges()
    
    for enlace in enlaces:
        
        nodo1, nodo2 = enlace
        
        if red_enfermita.get_edge_data(nodo1, nodo2)['weight'] < tol:
            
            red_enfermita.remove_edge(nodo1, nodo2)
            
    return red_enfermita

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
    
tolerancia = 15

solo_reacters = proyect_sin_posters(red_completa)

proyectada = proyect_enfermitos(solo_reacters, tolerancia)

nx.write_gml(proyectada, save_path_red+'tolerancia_{}.gml'.format(tolerancia)) # Red proyectada completa       


# Calculamos la componente gigante:
gc = proyectada.copy() 

gc_nodes = max(nx.connected_components(proyectada), key = len)

nodos_no_gc = set(gc.nodes()) - set(gc_nodes)

gc.remove_nodes_from(nodos_no_gc)

nx.write_gml(gc, save_path_red+'cg_{}.gml'.format(tolerancia))

#%%

tolerancia = 10

red_enf = nx.read_gml('save_path_red'+'cg_{}.gml')






















'''
vamos a proyectar primero a una red de reaccionadores, luego sobre esa red hacemos la proyeccion pesada y usamos el 
peso de cada enlace (cantidad de vecinos compartidos) como parametro a colocarle una tolerancia.

Despues podemos usar la proyeccion de la red dirigida sobre los enfermitos para encontrar (nuevamente mediante el peso)
la cantidad de enfermitos que likearon varios posts de un posteador.
'''

'''
NOTAS:

establecer un vector para cada enfermito que determine a cuantos post de cada categoria le dio like y otro para posteos    
proyectar sobre la red de enfermitos. 
ver comunidades por infomap, fast greedy, etc.
hacer dendograma por modularidad

2 analisis posibles: hacer silhouette de las comunidades vs las categorias que establecimos
o hacer un promedio d

hacer proyeccion incluyendo posters y ver si hay enfermitos que atraen muchos likes (enfermitos "famosos")


# Mirar pagina: https://networkx.github.io/documentation/stable/reference/algorithms/bipartite.html en la parte
# weightened projections y mirar la que habla del ratio de vecinos compartidos. Ver si eso sirve para establecer una
#tolerancia.

2) RED DE POPULARIDAD ENFERMITA:
    
    una proyeccion en la que se añade un enlace entre un posteador y un reacter si el reacter le dio like a n o mas posts
    suyos. (podria ser dirigida de posteador a reacter)
    
    Pintar a los enfermitos con la categoria que le toco en el analisis de comunidades y ver si hay homofilia (enfermitos
   de igual categoria tienen a estar juntos)
    
'''




