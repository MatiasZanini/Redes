# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 17:50:08 2020

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
import analisis_final_mati as func

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





















