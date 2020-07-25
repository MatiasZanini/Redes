# -*- coding: utf-8 -*-
"""
Created on Thu Jul  9 16:54:37 2020

@author: Mati
"""
# CARGA DE PAQUETES


#import os
#from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import networkx as nx
#from networkx.algorithms import bipartite
import igraph as ig
import matplotlib.cm as cm
import community as community_louvain
#import pickle
import funciones_final_mati as func
#import random

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

tolerancias = [2,3,4,6,7,8,9,11,12,13,14]

#tolerancias = [1]

red_completa.remove_nodes_from(tibios)

for tolerancia in tolerancias:
        
    
    
    #tolerancia = 12
    
    solo_reacters = func.proyect_sin_posters(red_completa)
    
    proyectada = func.proyect_enfermitos(solo_reacters, tolerancia)
    
    nx.write_gml(proyectada, save_path_red+'tolerancia_{}.gml'.format(tolerancia)) # Red proyectada completa
    
    #nx.write_gml(proyectada, save_path_red+'tolerancia_tibios_{}.gml'.format(tolerancia)) # Red proyectada completa       
    
    
    # Calculamos la componente gigante:
    gc = proyectada.copy() 
    
    gc_nodes = max(nx.connected_components(proyectada), key = len)
    
    nodos_no_gc = set(gc.nodes()) - set(gc_nodes)
    
    gc.remove_nodes_from(nodos_no_gc)
    
    nx.write_gml(gc, save_path_red+'cg_{}.gml'.format(tolerancia))
    
    #nx.write_gml(gc, save_path_red+'cg_tibios_{}.gml'.format(tolerancia))
    

#%% ---------------------------- CARGAR COMPONENTE GIGANTE ---------------------------------------------

tolerancia = 4 # Poner la tolerancia utilizada

red_enf = nx.read_gml(save_path_red+'cg_{}.gml'.format(tolerancia))

#red_enf = nx.read_gml(save_path_red+'cg_tibios_{}.gml'.format(tolerancia))

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

lou_dict = community_louvain.best_partition(red_enf, random_state = 10)


'''
Estimación de las modularidades:
'''

Q_infomap = red_enf_ig.modularity(infomap_ig)

Q_fg = red_enf_ig.modularity(fg_ig)

Q_lou = func.calcular_modularidad(red_enf, lou_dict)

print('La modularidad por infomap es:', Q_infomap)

print('La modularidad por Fast Greedy es:', Q_fg)

print('La modularidad por Louvain es:', Q_lou)

#%% ------------------------------- AGREGAMOS LA COMUNIDAD COMO ATRIBUTO A LOS NODOS ------------------

save_path_red = 'C:/Users/Mati/Documents/GitHub/Redes/datos/redes_gml/'

red = red_enf

nx.set_node_attributes(red, fg_dict, 'Comuna Fast Greedy')

nx.set_node_attributes(red, lou_dict, 'Comuna Louvain')

nx.write_gexf(red_enf, save_path_red + 'comunas_tol_{}.gexf'.format(tolerancia))

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
    
    comunidades.append([i for i in comu_dict if comu_dict[i]==comunidad]) # Ponemos las comunidades en cada elemento de la lista

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
    
plt.ylabel('Proporción de Likes', fontsize = 20)

plt.title('Distribución real de likes por comuna', fontsize = 20)

plt.legend(fontsize=16)

plt.xticks([1,2,3,4], labels = ['Categoría 1', 'Categoría 2', 'Categoría 3', 'Categoría 4'], fontsize=14)

plt.grid()


#%% ------------------------------------- COMPARACION CON MODELO NULO -------------------------------------------

'''
Respetando la distribucion de likes en cada categoria, asignarle una categoria preferida a cada enfermito y repetir
el analisis (modelo nulo). Si los histogramas de cada comunidad dan lo mismo, entonces descartamos la hipotesis. Sino,
hay algo.
'''

distrib = [] # Distribucion real total de categorias, dada por los datos.

for enf in red_enf.nodes():

    distrib.append(cat_norm[enf])

distrib = np.asarray(distrib)

distrib = sum(distrib)/sum(sum(distrib))

#%%

cat_norm_nulo = {} # Diccionario que asigna a cada enfermito una categoria al azar, respetando la distribucion al azar.

for enf in red_enf.nodes():
    
    cat_norm_nulo[enf] = int(func.dado_cargado([1,2,3,4], distrib))


#%%

comunidades_dict_nulo = [] # Lista de diccionarios enf:categoria para cada comuna

categorias = [1,2,3,4]

for comunidad in comunidades_numero:
    
    comunidades_dict_nulo.append({key : cat_norm_nulo[key] for key in comunidades[comunidad]})


cat_por_com_nulo = [] # Lista con la distribucion nula para cada comunidad

for comunidad in comunidades_dict_nulo:
    
    counts, bins = np.histogram(list(comunidad.values()), bins = categorias+[len(categorias)+1])

    cat_por_com_nulo.append(counts/sum(counts))
    
    
#%% -------------------------- MODELO NULO PROMEDIADO --------------------------------------------

n_prom = 100

distribuciones_prom = []

for i in range(n_prom):
    
    cat_norm_nulo = {}
    
    for enf in red_enf.nodes():
        
        cat_norm_nulo[enf] = int(func.dado_cargado([1,2,3,4], distrib))

    comunidades_dict_nulo = [] # Lista de diccionarios enf:categoria para cada comuna
    
    categorias = [1,2,3,4]
    
    for comunidad in comunidades_numero:
        
        comunidades_dict_nulo.append({key : cat_norm_nulo[key] for key in comunidades[comunidad]})
    
    
    cat_por_com_nulo = [] # Lista con la distribucion nula para cada comunidad
    
    for comunidad in comunidades_dict_nulo:
        
        counts, bins = np.histogram(list(comunidad.values()), bins = categorias+[len(categorias)+1])
    
        cat_por_com_nulo.append(counts/sum(counts))
        
        
    distribuciones_prom.append(np.asarray(cat_por_com_nulo))


cat_por_com_nulo_prom = []

for fila in sum(distribuciones_prom):
    
    cat_por_com_nulo_prom.append(fila/sum(fila))    
    


#%% Grafico

cont = 0

for i in cat_por_com_nulo_prom:
    
    plt.plot([1,2,3,4], i, '.-', label= 'Comunidad {}'.format(cont))
    
    cont += 1
    
#plt.xlabel('Categorias', fontsize = 20)

plt.ylabel('Proporción de likes', fontsize = 20)

plt.title('Distribución de likes por comuna en el Modelo Nulo 1', fontsize = 20)

plt.legend(fontsize=16)

plt.xticks([1,2,3,4], labels = ['Categoría 1', 'Categoría 2', 'Categoría 3', 'Categoría 4'], fontsize=14)

plt.grid()


#%% ---------------------------- MODELO NULO 2 ------------------------------------------------------------

'''
La idea es crear la red mezclando las categorias de los post de forma aleatoria.
'''

props_nulo_2 = {}

cat_norm_nulo_2 = {}

categorias = [1,2,3,4] #(humor negro, verde, serie, interno)

iter_enfermitos = enfermitos # Poner la particion en comunidades como diccionario o lista de enfermitos

for enf in iter_enfermitos:
    
    props_nulo_2[enf] = func.prop_enfermitos_nulo(enf, df_completo) # Diccionario con vectores de categorias por enfermito

for enf in props_nulo_2:
    
    cat_norm_nulo_2[enf] = func.categorias_norm(props_nulo_2[enf], categorias)
    
    
comunidades_dict_nulo_2 = [] # Lista de diccionarios enf:categoria para cada comuna

categorias = [1,2,3,4]

for comunidad in comunidades_numero:
    
    comunidades_dict_nulo_2.append({key : cat_norm_nulo_2[key] for key in comunidades[comunidad]})

#%% ------------------------------ GUARDAR MODELO NULO 2 -----------------------------------

nombre_algo = 'fg'

name_nulo = 'modelo_nulo_{}_tol_{}'.format(nombre_algo, tolerancia)

func.save_dict(comunidades_dict_nulo_2, path_props, name_nulo)

    
#%% ------------------------------- CARGAR MODELO NULO 2 -------------------------------------

'''
Si ya se realizó el modelo nulo para una dada red, cargarlo desde aca.

NOTA: lo correcto sería hacerlo de cero cada vez.
'''

nombre_algo = 'fg'

name_nulo = 'modelo_nulo_{}_tol_{}'.format(nombre_algo, tolerancia)


comunidades_dict_nulo_2 = func.load_dict(path_props, name_nulo)



#%%

cat_por_com_nulo_2 = []

for comunidad2 in comunidades_dict_nulo_2:
    
    cat_por_com_nulo_2.append(sum(comunidad2.values())/sum(sum(comunidad2.values()))) # Lista con los pesos de cada categoria en cada comuna

cont = 0

for i in cat_por_com_nulo_2:
    
    plt.plot([1,2,3,4], i, '.-', label= 'Comunidad {}'.format(cont))
    
    cont += 1
    
plt.ylabel('Proporción de likes', fontsize = 20)

plt.title('Distribución de likes por comuna en el Modelo Nulo 2', fontsize = 20)

plt.legend(fontsize=16)

plt.xticks([1,2,3,4], labels = ['Categoría 1', 'Categoría 2', 'Categoría 3', 'Categoría 4'], fontsize=14)

plt.grid()




#%% ------------------------------ CORRELACION ENTRE DISTRIBUCION REAL Y MODELOS NULOS ------------------

'''
Aca pondremos todos los vectores juntos, calcularemos la correlacion-->similitud--> distancia y haremos un dendograma 
comparandolos.
'''
cmap = cm.get_cmap('RdYlGn')

labels = ['Com. 0', 'Com. 1', 'Com. 2', 'Com. 3', 'Com. 0 M.N.', 'Com. 1 M.N.', 
          'Com. 2 M.N.', 'Com. 3 M.N.']


# MODELO NULO 1

corr1 = np.corrcoef(cat_por_com + cat_por_com_nulo)

sim1 = (1 + corr1 )/2


fig1, ax1 = plt.subplots(figsize=(20,20))

cax1 = ax1.matshow(sim1, cmap = cmap)

plt.title('Similaridad entre comunas y Modelo Nulo 1', fontsize=20)

plt.xticks(range(8), labels)

plt.yticks(range(8), labels)

fig1.colorbar(cax1, ticks=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, .75,.8,.85,.90,.95,1])

plt.show()



#%%
# MODELO NULO 2

corr2 = np.corrcoef(cat_por_com + cat_por_com_nulo_2)

sim2 = (1 + corr2 )/2


fig2, ax2 = plt.subplots(figsize=(20,20))

cax2 = ax2.matshow(sim2, cmap = cmap)

plt.title('Similaridad entre comunas y Modelo Nulo 2', fontsize=20)

plt.xticks(range(8), labels)

plt.yticks(range(8), labels)

fig2.colorbar(cax2, ticks=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, .75,.8,.85,.90,.95,1])

plt.show()


#%%

# SIMILITUD ENTRE ENFERMITOS DE UNA MISMA COMUNIDAD

enf_ordenados = []

for comuna in comunidades:
    
    enf_ordenados = enf_ordenados + comuna # Lista con los enfermitos ordenados por comuna

cat_ordenadas = []

for enf in enf_ordenados:
    
    cat_ordenadas.append(cat_norm[enf]) # Lista con los vectores de categorias likeadas por cada enfermito

corr_enf = np.corrcoef(cat_ordenadas)

sim_enf = (1+corr_enf)/2


ticks = []

tick_count = 0

for i in range(len(comunidades)):
    
    ticks.append(tick_count)
    
    tick_count += len(comunidades[i])

label_com = []

for i in range(len(comunidades)):
    
    label_com.append('Inicio Comunidad {}'.format(i))



fig_enf, ax_enf = plt.subplots(figsize=(20,20))

cax_enf = ax_enf.matshow(sim_enf, cmap = cmap)

plt.title('Similaridad entre Reaccionadores', fontsize=20)

plt.xticks(ticks, label_com)

plt.yticks(ticks, label_com)

fig_enf.colorbar(cax_enf, ticks=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, .75,.8,.85,.90,.95,1])

plt.show()





# labels = []
# for hood in hood_menu_data:
# labels.append(hood["properties"]['NAME'])
 
# fig, ax = plt.subplots(figsize=(20,20))
# cax = ax.matshow(hood_cosine_matrix, interpolation='nearest')
# ax.grid(True)
# plt.title('San Francisco Similarity matrix')
# plt.xticks(range(33), labels, rotation=90);
# plt.yticks(range(33), labels);
# fig.colorbar(cax, ticks=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, .75,.8,.85,.90,.95,1])
# plt.show()





#%% ------------------ HISTOGRAMAS PARA CATEGORIAS PREFERIDAS (analisis viejo) -------------------------


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




