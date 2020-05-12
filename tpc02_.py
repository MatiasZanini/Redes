# -*- coding: utf-8 -*-
"""
                **Introducción a Redes Complejas en Biología de Sistemas**
                        Trabajo Computacional 2 - Estructura a Gran Escala 
Entrega 14/05


Grupo: Camila Sanz, Matías Zanini.
"""

################################################################################
#                                 PAQUETES 
################################################################################

#import pandas as pd       #DB
import numpy as np        #matemática, simil maltab 
import networkx as nx
import matplotlib.pyplot as plt
#from matplotlib_venn import venn3
#import igraph as ig
import random
import math
# Evitar acumulación de mensajes de warning en el display
import warnings  
warnings.filterwarnings("ignore")


#%%

################################################################################
#                               PUNTO 1 
################################################################################
















#%%

################################################################################
#                               PUNTO 2 
################################################################################

#Definimos las funciones que vamos a usar a lo largo del punto:
    
def moneda_cargada(p):
    
    '''
    Devuelve un booleano que indica True si 
    '''
    
    cara = random.random()

    if cara <= p:
        
        return True
    
    else:
        
        return False


'''
a)
'''

# Armamos la red random cuyos nodos se conectan con probabilidad constante p:
    
p = 0.2 # Probabilidad con la que se conectan dos nodos.

n = 1e4 # Cantidad de Nodos en la red

nodos = np.arange(1, n+1, 1)

red_er = nx.Graph() # Inicializamos el grafo

red_er.add_nodes_from(nodos) 

enlaces_list = [] 

for fila in nodos:
    
    #Nos movemos por arriba de la diagonal de la matriz de adyacencia para no contar 2 veces los pares de nodos:
    for columna in range(int(fila)+1, len(nodos)+1): 
    
        if moneda_cargada(p):
            
            # Añadimos a la lista de enlaces la tupla si sale favorable la tirada con probabilidad p:
            enlaces_list.append( (fila, columna) ) 
            
red_er.add_edges_from(enlaces_list)

    
    























#%%

################################################################################
#                               PUNTO 3 
################################################################################