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
import funciones_final_mati as func


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
#                                      ANALISIS
#####################################################################################


# ------------------------------- PROPIEDADES DE LOS ENFERMITOS -----------------------------

props = {}

tag = {}

cat_norm = {}

categorias = [1,2,3,4] #(humor negro, verde, serie, interno)

iter_enfermitos = enfermitos # Poner la particion en comunidades como diccionario o lista de enfermitos

for enf in iter_enfermitos:
    
    props[enf] = func.prop_enfermitos(enf, df_completo) # Diccionario con vectores de categorias por enfermito

for enf in props:
    
    tag[enf] = func.tag_enfermitos(enf, df_completo, categorias) # Diccionario con categoria favorita por enfermito

for enf in props:
    
    cat_norm[enf] = func.categorias_norm(props[enf], categorias)
        

#%% ------------------ GUARDAR LAS PROPIEDADES DE LOS ENFERMITOS ------------------------------

path_props = 'C:/Users/Mati/Documents/GitHub/Redes/datos/props_enfermitos/'

name_tag = 'tags_enfermitos_sin_cat_5'

name_props = 'propiedades_enfermitos_sin_cat_5'

name_cat_norm = 'cat_normalizadas_sin_cat_5'

func.save_dict(props, path_props, name_props)

func.save_dict(tag, path_props, name_tag)

func.save_dict(cat_norm, path_props, name_cat_norm)


#%% -------------------- CARGAR LAS PROPIEDADES DE LOS ENFERMITOS -----------------------------

path_props = 'C:/Users/Mati/Documents/GitHub/Redes/datos/props_enfermitos/'

name_tag = 'tags_enfermitos_sin_cat_5'

name_props = 'propiedades_enfermitos_sin_cat_5'

name_cat_norm = 'cat_normalizadas_sin_cat_5'

props = func.load_dict(path_props, name_props)

tag = func.load_dict(path_props, name_tag)

cat_norm = func.load_dict(path_props, name_cat_norm)





#%%----------------------------------- PROYECCIONES ----------------------------------

tibios = []

for enf in tag:
    
    tipo_array = type(np.array([]))
    
    if type(tag[enf]) == tipo_array:
        
        tibios.append(enf)
        
#%% ---------------------------- GUARDAR RED ANALIZADA ------------------------------

#tolerancias = [2,3,4,6,7,8,9,11,12,13,14]

tolerancias = [1]

red_completa.remove_nodes_from(tibios)

for tolerancia in tolerancias:
        
    
    
    #tolerancia = 12
    
    solo_reacters = func.proyect_sin_posters(red_completa)
    
    proyectada = func.proyect_enfermitos(solo_reacters, tolerancia)
    
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

pesos = {}

for enlace in red_enf.edges():
    
    u,v = enlace
    
    pesos[enlace]=red_enf.get_edge_data(u,v)['weight']

red_enf_ig.es['weight'] = list(pesos.values())


# Infomap:

infomap_ig = red_enf_ig.community_infomap(edge_weights='weight') # Creamos la partición mediante el algoritmo de Infomap

infomap_dict = func.ig_part2dict(red_enf_ig, infomap_ig) # Convertimos la partición por Infomap a un diccionario


# Fast Greedy:

fg_dendograma = red_enf_ig.community_fastgreedy(weights= 'weight') # Dendograma con el FastGreedy que maximiza Q

fg_ig = fg_dendograma.as_clustering() # Comunidades dadas por Fast Greedy como objeto de Igraph

fg_dict = func.ig_part2dict(red_enf_ig, fg_ig)

lou_dict = community_louvain.best_partition(red_enf)


'''
Estimación de las modularidades:
'''

Q_infomap = red_enf_ig.modularity(infomap_ig)

Q_fg = red_enf_ig.modularity(fg_ig)

Q_lou = func.calcular_modularidad(red_enf, lou_dict)

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

func.graficar_particion(red_enf, infomap_dict, posiciones,label='Categoría Preferida', tamano_nodo=tamano_nodo, colormap=colormap)

plt.subplot(1,3,2)

plt.title('Comunidades por Fast Greedy')

func.graficar_particion(red_enf, fg_dict, posiciones,label='Categoría Preferida', tamano_nodo=tamano_nodo, colormap=colormap)

plt.subplot(1,3,3)

plt.title('Comunidades por Louvain')

func.graficar_particion(red_enf, lou_dict, posiciones, label='Categoría Preferida', tamano_nodo = tamano_nodo, colormap =colormap)



#%% -------------------------------- CLASIFICANDO COMUNIDADES ---------------------------------------------

comu_dict = fg_dict # Poner el diccionario con las comunidades

comunidades_numero = list(set(comu_dict.values()))

comunidades = []

for comunidad in comunidades_numero:
    
    comunidades.append([i for i in lou_dict if lou_dict[i]==comunidad]) # Ponemos las comunidades en cada elemento de la lista

comunidades_dict = []

comunidades_tag_dict = []

for comunidad in comunidades_numero:
    
    comunidades_tag_dict.append({key : tag[key] for key in comunidades[comunidad]}) # Lista con diccionario de las comus.
    
for comunidad in comunidades_numero:
    
    comunidades_dict.append({key : cat_norm[key] for key in comunidades[comunidad]}) # Lista con diccionario de las comus.


#%%

cat_por_com = []

for comunidad in comunidades_dict:
    
    cat_por_com.append(sum(comunidad.values())/sum(sum(comunidad.values()))) # Lista con los pesos de cada categoria en cada comuna
cont = 0
for i in cat_por_com:
    
    plt.plot([1,2,3,4], i, '.-', label= 'Comunidad {}'.format(cont))
    
    cont += 1
    
plt.xlabel('Categorias')

plt.ylabel('Cantidad de likes')

plt.legend()


#%% ------------------------------------- COMPARACION CON MODELO NULO -------------------------------------------

'''
Respetando la distribucion de likes en cada categoria, asignarle una categoria preferida a cada enfermito y repetir
el analisis (modelo nulo). Si los histogramas de cada comunidad dan lo mismo, entonces descartamos la hipotesis. Sino,
hay algo.
'''


































#%% ------------------ HISTOGRAMAS PARA CATEGORIAS PREFERIDAS -------------------------


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




