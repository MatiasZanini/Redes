# -*- coding: utf-8 -*-
"""
Created on Thu Jul  9 16:54:37 2020

@author: Mati
"""
# CARGA DE PAQUETES


import os
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import networkx as nx
from networkx.algorithms import bipartite
import igraph as ig
import matplotlib.cm as cm
import community as community_louvain
import pickle

#%%


#####################################################################################
#                               UNIFICACION DE DATAFRAMES
#####################################################################################


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

####################################################################################
#                               SIN CATEGORIA 5
####################################################################################

# Categorías: [humor negro, humor verde, humor de serie, humor interno] (del 1 al 4)


save_path = 'C:/Users/Mati/Documents/GitHub/Redes/datos/data_frames/'

df_total_name = 'data_total.p'

df_parcial_name = 'data_sin_cat_5.p'

data_parcial = pd.read_pickle(save_path + df_total_name)

df1 = pd.read_pickle(save_path + 'primera_tanda.p')

df2 = pd.read_pickle(save_path + 'segunda_tanda.p')

df3 = pd.read_pickle(save_path + 'datos_red.p')

borrables = [3,8,9,20,47,46, 48, 49, 50, 52, 54,56,53]

df3.drop(borrables, inplace = True) # Eliminamos posts mal categorizados.

df1.replace({'5' : '3', '6': '4', '7':'5'}, inplace=True)

for i in range(len(df1)):
    
    if df1['categoria'][i] == '5':
        
        df1.drop(index=i, inplace = True)

df2.replace({'5' : '3', '6': '4', '7':'5'}, inplace=True)

data_parcial = df1

data_parcial = data_parcial.append(df2, ignore_index=True)

data_parcial = data_parcial.append(df3, ignore_index=True)

data_parcial.drop_duplicates(subset ='url', keep = 'first', inplace = True)

data_parcial.to_pickle(save_path + df_parcial_name)


#%%
#####################################################################################
#                               CARGA DE ENFERMITOS Y POSTS
#####################################################################################

save_path = 'C:/Users/Mati/Documents/GitHub/Redes/datos/data_frames/'

save_path_red = 'C:/Users/Mati/Documents/GitHub/Redes/datos/redes_gml/'

#save_name = 'data_total.p'

save_name = 'data_sin_cat_5.p'

#name_red = 'red_completa.gml'

name_red = 'red_sin_cat_5.gml'

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


#%% ------------------------------- GENERAMOS LA RED Y LA GUARDAMOS -------------------------------

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


nx.write_gml(red_completa, save_path_red + name_red)

#%% 

#####################################################################################
#                                CARGA DE RED COMPLETA Y DATAFRAME
#####################################################################################

# Carga de datos

save_path = 'C:/Users/Mati/Documents/GitHub/Redes/datos/data_frames/'

save_path_red = 'C:/Users/Mati/Documents/GitHub/Redes/datos/redes_gml/'

#red_name = 'red_completa.gml'

red_name = 'red_sin_cat_5.gml'

#df_name = 'data_total.p'

df_name = 'data_sin_cat_5.p'

red_completa = nx.read_gml(save_path_red + red_name)

df_completo = pd.read_pickle(save_path + df_name)

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



          
#%%

#####################################################################################
#                                      ANALISIS
#####################################################################################


# ------------------------------- PROPIEDADES DE LOS ENFERMITOS -----------------------------

props = {}

tag = {}

categorias = [1,2,3,4] #(humor negro, verde, serie, interno)

iter_enfermitos = enfermitos # Poner la particion en comunidades como diccionario o lista de enfermitos

for enf in iter_enfermitos:
    
    props[enf] = prop_enfermitos(enf, df_completo) # Diccionario con vectores de categorias por enfermito

for enf in props:
    
    tag[enf] = tag_enfermitos(enf, df_completo, categorias)
        
# tag_df = pd.DataFrame(index = range(len(iter_enfermitos)),columns = ['enfermito', 'categoria pref'])

# for enf in range(len(iter_enfermitos)):
    
#     tag_df['categoria pref'][enf] = tag[enfermitos[enf]]

#     tag_df['enfermito'][enf] = enfermitos[enf]

#%% ------------------ GUARDAR LAS PROPIEDADES DE LOS ENFERMITOS ------------------------------

path_props = 'C:/Users/Mati/Documents/GitHub/Redes/datos/props_enfermitos/'

name_tag = 'tags_enfermitos_sin_cat_5'

name_props = 'propiedades_enfermitos_sin_cat_5'

save_dict(props, path_props, name_props)

save_dict(tag, path_props, name_tag)


#%% -------------------- CARGAR LAS PROPIEDADES DE LOS ENFERMITOS -----------------------------

path_props = 'C:/Users/Mati/Documents/GitHub/Redes/datos/props_enfermitos/'

name_tag = 'tags_enfermitos_sin_cat_5'

name_props = 'propiedades_enfermitos_sin_cat_5'

props = load_dict(path_props, name_props)

tag = load_dict(path_props, name_tag)



#%%----------------------------------- PROYECCIONES ----------------------------------

tibios = []

for enf in tag:
    
    tipo_array = type(np.array([]))
    
    if type(tag[enf]) == tipo_array:
        
        tibios.append(enf)
        
#%% ---------------------------- GUARDAR RED ANALIZADA ------------------------------

tolerancias = [2,3,4,6,7,8,9,11,12,13,14]

red_completa.remove_nodes_from(tibios)

for tolerancia in tolerancias:
        
    
    
    #tolerancia = 12
    
    solo_reacters = proyect_sin_posters(red_completa)
    
    proyectada = proyect_enfermitos(solo_reacters, tolerancia)
    
    nx.write_gml(proyectada, save_path_red+'tolerancia_{}.gml'.format(tolerancia)) # Red proyectada completa       
    
    
    # Calculamos la componente gigante:
    gc = proyectada.copy() 
    
    gc_nodes = max(nx.connected_components(proyectada), key = len)
    
    nodos_no_gc = set(gc.nodes()) - set(gc_nodes)
    
    gc.remove_nodes_from(nodos_no_gc)
    
    nx.write_gml(gc, save_path_red+'cg_{}.gml'.format(tolerancia))
    

#%% ---------------------------- CARGAR COMPONENTE GIGANTE ---------------------------------------------

tolerancia = 11 # Poner la tolerancia utilizada

red_enf = nx.read_gml(save_path_red+'cg_{}.gml'.format(tolerancia))

nx.set_node_attributes(red_enf, tag, 'Categoría Preferida')


red_enf_ig = ig.Graph.TupleList(red_enf.edges(), directed=False)


#%% ---------------------------- AGREGARLE EL TAG A LOS NODOS COMO PROPIEDAD ---------------------------

#red = red_completa

red = red_enf

nx.set_node_attributes(red, tag, 'Categoría Preferida')





#%%

#---------------------------------- CALCULO DE COMUNIDADES -------------------------------------------

# Infomap:

infomap_ig = red_enf_ig.community_infomap() # Creamos la partición mediante el algoritmo de Infomap

infomap_dict = ig_part2dict(red_enf_ig, infomap_ig) # Convertimos la partición por Infomap a un diccionario


# Fast Greedy:

fg_dendograma = red_enf_ig.community_fastgreedy(weights=None) # Dendograma con el FastGreedy que maximiza Q

fg_ig = fg_dendograma.as_clustering() # Comunidades dadas por Fast Greedy como objeto de Igraph

fg_dict = ig_part2dict(red_enf_ig, fg_ig)

lou_dict = community_louvain.best_partition(red_enf)


'''
Estimación de las modularidades:
'''

Q_infomap = red_enf_ig.modularity(infomap_ig)

Q_fg = red_enf_ig.modularity(fg_ig)

Q_lou = calcular_modularidad(red_enf, lou_dict)

print('La modularidad por infomap es:', Q_infomap)

print('La modularidad por Fast Greedy es:', Q_fg)

print('La modularidad por Louvain es:', Q_lou)


#%%

#---------------------------------- GRAFICO DE PARTICIONES --------------------

posiciones = nx.kamada_kawai_layout(red_enf)

colormap = 'Accent'

tamano_nodo = 35

plt.figure(figsize=(15,15))

plt.subplot(1,3,1)

plt.title('Comunidades por Infomap')

graficar_particion(red_enf, infomap_dict, posiciones,label='Categoría Preferida', tamano_nodo=tamano_nodo, colormap=colormap)

plt.subplot(1,3,2)

plt.title('Comunidades por Fast Greedy')

graficar_particion(red_enf, fg_dict, posiciones,label='Categoría Preferida', tamano_nodo=tamano_nodo, colormap=colormap)

plt.subplot(1,3,3)

plt.title('Comunidades por Louvain')

graficar_particion(red_enf, lou_dict, posiciones, label='Categoría Preferida', tamano_nodo = tamano_nodo, colormap =colormap)



#%% -------------------------------- CLASIFICANDO COMUNIDADES ---------------------------------------------

comu_dict = lou_dict # Poner el diccionario con las comunidades

comunidades_numero = list(set(comu_dict.values()))

comunidades = []

for comunidad in comunidades_numero:
    
    comunidades.append([i for i in lou_dict if lou_dict[i]==comunidad]) # Ponemos las comunidades en cada elemento de la lista

comunidades_dict = []

for comunidad in comunidades_numero:
    
    comunidades_dict.append({key : tag[key] for key in comunidades[comunidad]}) # Lista con diccionario de las comus.

#%%


comunidad0 = list(comunidades_dict[0].values())

bins = [1,2,3,4,5]

counts0, bins0 = np.histogram(comunidad0, bins=bins)

comunidad1 = list(comunidades_dict[1].values())

bins = [1,2,3,4,5]

counts1, bins1 = np.histogram(comunidad1, bins=bins)

comunidad2 = list(comunidades_dict[2].values())

bins = [1,2,3,4,5]

counts2, bins2 = np.histogram(comunidad2, bins=bins)

comunidad3 = list(comunidades_dict[3].values())

bins = [1,2,3,4,5]

counts3, bins3 = np.histogram(comunidad3, bins=bins)

plt.figure() # Sin normalizar

plt.plot([1,2,3,4], counts0, '.-', label= 'Comunidad 0')

plt.plot([1,2,3,4], counts1, '.-', label= 'Comunidad 1')

plt.plot([1,2,3,4],counts2, '.-', label= 'Comunidad 2')

plt.plot([1,2,3,4],counts3, '.-', label= 'Comunidad 3')

plt.xlabel('Categorias')

plt.ylabel('Cantidad de likes')

plt.legend()

#%%

plt.figure() # Normalizado

plt.plot([1,2,3,4], counts0/max(counts0), '.-', label= 'Comunidad 0')

plt.plot([1,2,3,4], counts1/max(counts1), '.-', label= 'Comunidad 1')

plt.plot([1,2,3,4],counts2/max(counts2), '.-', label= 'Comunidad 2')

plt.plot([1,2,3,4],counts3/max(counts3), '.-', label= 'Comunidad 3')

plt.xlabel('Categorias')

plt.ylabel('Cantidad de enfermitos en la categoría')

plt.legend()




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




