# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 17:35:30 2020

@author: Mati
"""
import pickle
import numpy as np
import networkx as nx
from networkx.algorithms import bipartite
import matplotlib.cm as cm

#%%
#####################################################################################
#                                    FUNCIONES
#####################################################################################

def save_dict(dictionary, path, name):
    
    with open(path + name + '.pkl', 'wb') as f:
        
        pickle.dump(dictionary, f, pickle.HIGHEST_PROTOCOL)

def load_dict(path, name):
    
    with open(path + name + '.pkl', 'rb') as f:
        
        return pickle.load(f)


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
    si 2 enfermitos le dieron like a "tol" o mas posts iguales --> enlace pesado por la cantidad de posts en comun
    '''
    
    enfermitos, posts = nx.bipartite.sets(red)
    
    red_enfermita = bipartite.weighted_projected_graph(red, enfermitos)
    
    enlaces = red_enfermita.edges()
    
    for enlace in enlaces:
        
        nodo1, nodo2 = enlace
        
        if red_enfermita.get_edge_data(nodo1, nodo2)['weight'] < tol:
            
            red_enfermita.remove_edge(nodo1, nodo2)
            
    return red_enfermita



def proyect_enfermitos_sin_peso(red, tol):
    
    '''
    si 2 enfermitos le dieron like a "tol" o mas posts iguales --> enlace sin peso
    '''
    
    enfermitos, posts = nx.bipartite.sets(red)
    
    red_enfermita_w = bipartite.weighted_projected_graph(red, enfermitos)
    
    red_enfermita = nx.Graph()
    
    enlaces = red_enfermita_w.edges()
    
    for enlace in enlaces:
        
        nodo1, nodo2 = enlace
        
        if red_enfermita.get_edge_data(nodo1, nodo2)['weight'] >= tol:
            
            red_enfermita.add_edge(nodo1, nodo2)
            
    return red_enfermita



def graficar_particion(Red, particion_diccionario, posiciones, label=None, tamano_nodo = 200, colormap = 'viridis'):
    
    '''
    Grafica una red de networkx con colores segun su partició dada por un diccionario según un dado set de posiciones 
    para los nodos.
    '''
    
    if label:
        
        labels = nx.get_node_attributes(Red, label)
    
    
    cmap = cm.get_cmap(colormap, max(particion_diccionario.values())+1) # viridis es el mapa de colores
    
    # grafico los nodos
    nx.draw_networkx_nodes(Red, posiciones, particion_diccionario.keys(), labels = labels, node_size = tamano_nodo,
                           cmap=cmap, node_color = list(particion_diccionario.values()), with_labels = False)
    
    # grafico los enlaces aparte
    nx.draw_networkx_edges(Red, pos = posiciones, alpha=0.5)
    
    

def ig_part2dict(Red_igraph, particion_igraph):
    
    '''
    Convierte una partición en comunidades a un diccionario
    '''
    
    particion_dict = {}
    
    for cluster in range(len(particion_igraph)):
    
        for nodo in Red_igraph.vs(particion_igraph[cluster])['name']:
        
            particion_dict.update({nodo:cluster})
  
    return particion_dict

def prop_enfermitos(enf, df):
    
    '''
    Pide un enfermito y el dataframe con los datos de posts, posteadores y reacters. 
    Devuelve una lista con las categorias de los posts que likeo
    '''
    
    prop = []
    
    for i in range(len(df)):
        
        if enf in df['reacters'][i]:
            
            prop.append(int(df['categoria'][i]))
    
    
    return prop


def tag_enfermitos(enf, df, categorias):
    
    '''
    Pide un enfermito, el dataframe con los datos de los posts y una lista con las categorias.
    Devuelve un int con la categoria preferida del enfermito o un array con las categorias preferidas si era mas de una.
    '''
    
    
    counts, bins = np.histogram(np.asarray(prop_enfermitos(enf, df)), bins = categorias+[len(categorias)+1])
    
    index_maximos = [i for i, x in enumerate(counts) if x == max(counts)]
    
    if len(index_maximos)==1:
        
        tag = bins[index_maximos[0]]
        
    else:
        
        tag = bins[index_maximos] # Si tiene mas de una categoria preferida, se guarda un array con dichas categorias.

    return tag


def categorias_norm(likes, categorias):
    
    cant_posts = len(likes)
    
    counts, bins = np.histogram(likes, bins = categorias+[len(categorias)+1])
    
    cat_norm = np.asarray(counts/cant_posts)
    
    return cat_norm




def armar_listas_por_comunidades(particion_diccionario):
  
    listas_nodos_particiones = []
  
    comunidades = 0
  
    lista_comunidad = []
  
    for nodo in particion_diccionario.keys():
    
        if particion_diccionario[nodo] == comunidades:
      
            lista_comunidad.append(nodo)
    
        else:
        
            comunidades+=1
          
            listas_nodos_particiones.append(lista_comunidad)
          
            lista_comunidad = [nodo]
  
    listas_nodos_particiones.append(lista_comunidad)
  
    return listas_nodos_particiones

def calcular_modularidad(Red,particion_diccionario):
  
    listas_por_comunidades = armar_listas_por_comunidades(particion_diccionario)
  
    enlaces = len(Red.edges())
  
    enlaces_en_comunidad = []
  
    suma_grados_comunidad_2 = []
  
    for i in range(len(listas_por_comunidades)):
    
        enlaces_en_comunidad.append(len(Red.subgraph(listas_por_comunidades[i]).edges()))
    
        suma_grados_comunidad_2.append(sum([grado for nodo, grado in Red.degree(listas_por_comunidades[i])])**2)
  
    return sum(enlaces_en_comunidad)/enlaces - (sum(suma_grados_comunidad_2)/(4*enlaces**2))