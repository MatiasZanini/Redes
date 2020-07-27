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
from matplotlib import pyplot as plt
from random import choices
import random
from sklearn.metrics import normalized_mutual_info_score as mis # Función para calcular información mutua
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


def prop_enfermitos_nulo(enf, df):
    
    '''
    Pide un enfermito y el dataframe con los datos de posts, posteadores y reacters. 
    Mezcla las categorias de los post y devuelve una lista con las categorias nuevas de los posts que likeo
    '''
    
    prop = []
    
    categorias = list(df['categoria'])
    
    random.shuffle(categorias)
    
    random.shuffle(categorias)
    
    random.shuffle(categorias)
    
    for i in range(len(df)):
        
        if enf in df['reacters'][i]:
            
            prop.append(int(categorias[i]))
    
    
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



def tag_enfermitos_nulo(enf, df, categorias, props):
    
    '''
    Pide un enfermito, el dataframe con los datos de los posts y una lista con las categorias.
    Devuelve un int con la categoria preferida del enfermito o un array con las categorias preferidas si era mas de una.
    '''
    
    
    counts, bins = np.histogram(np.asarray(props[enf]), bins = categorias+[len(categorias)+1])
    
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


def dado_cargado(poblacion, distribucion, cant_muestras = 1):
    
    return np.asarray(choices(poblacion, distribucion, k = cant_muestras))
    

def modularidad(particion, matrix_pesos):
    
    W = np.sum(matrix_pesos)
    
    fuerzas = [np.sum(matrix_pesos[n,:]) for n in np.arange(0, matrix_pesos.shape[0])]
    
    Q = 0
    
    for i in np.arange(np.min(particion),np.max(particion)+1):
    
        indices = np.where(np.array(particion)==i)[0]
        
        for j in indices:
        
            for k in indices:
            
                Q = Q + matrix_pesos[j][k] - fuerzas[j] * fuerzas[k] /(2 * W)
    
    Q=Q / (2 * W)
    
    return Q


def silhouette(Red,particion): # La función nos pide la Red y la partición como diccionario
  S=[]
  limites=[0]
  # Recorro los clusters sin repetir
  for cluster in set(particion.values()): 
    #Filtro los nodos que pertenecen a este cluster
    nodos_en_cluster = [nodo for (nodo, value) in particion.items() if value == cluster] 
    S_cluster=[]
    # Recorro los nodos del cluster
    for nodo in nodos_en_cluster:
      distancias_dentro=[]
      distancias_fuera=[]
      # Recorro los nodos del mismo cluster
      for nodo_en_cluster in nodos_en_cluster:
        if nodo != nodo_en_cluster:
          # Calculo y guardo la distancia, si no es consigo mismo
          distancias_dentro.append(nx.shortest_path_length(Red, source=nodo, target=nodo_en_cluster, weight=None))
      # Recorro los nodos de los otros clusters
      for nodo_fuera in Red.nodes():
        if particion[nodo_fuera] != cluster:
          # Calculo y guardo la distancia
          distancias_fuera.append(nx.shortest_path_length(Red, source=nodo, target=nodo_fuera, weight=None)) 
      # Calculo la distancia media para los del mismo cluster
      distancia_media_dentro=np.mean(distancias_dentro)
      # Calculo la distancia mínima para los nodos fuera del cluster
      distancia_min_fuera=np.min(distancias_fuera)
      # Calculo y guardo la Silhouette del nodo
      S_cluster.append((distancia_min_fuera-distancia_media_dentro)/np.max([distancia_min_fuera,distancia_media_dentro]))
    # Ordeno las Silhouette del mismo cluster por valor, para graficar lindo
    S_cluster=sorted(S_cluster)
    # Me guardo en qué nodo termina cada cluster, para graficar clusters por colores
    limites.append(len(S_cluster)+limites[-1])
    # Agrego las Silhouette de este cluster a la lista de todas
    S = S + S_cluster
  # Calculo la Silhouette media
  S_media = np.mean(S)
  # Grafico todas con colores por clusters
  for i in range(len(limites)-1):
    plt.plot(range(limites[i],limites[i+1]),S[limites[i]:limites[i+1]])
  plt.plot(range(0,limites[-1]),S_media*np.ones(limites[-1]))
  return S,S_media



def info_mutua(particion_1,particion_2, red):
    
    info_1 = []
    
    info_2 = []
    
    for nodo in red.nodes():
      
        info_1.append(particion_1[nodo])
      
        info_2.append(particion_2[nodo])
    
    info_1 = np.array(info_1)
    
    info_2 = np.array(info_2)
    
    return mis(info_1, info_2, average_method='arithmetic') #Esta es la función que calcula info mutua
    
    
    
    




