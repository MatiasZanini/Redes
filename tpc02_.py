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

#%%

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
    
    # Nos movemos por arriba de la diagonal de la matriz de adyacencia para no contar 2 veces los pares de nodos:
    for columna in range(int(fila)+1, len(nodos)+1): 
    
        if moneda_cargada(p):
            
            # Añadimos a la lista de enlaces la tupla si sale favorable la tirada con probabilidad p:
            enlaces_list.append( (fila, columna) ) 
            
red_er.add_edges_from(enlaces_list)

# i.

m =  red_er.number_of_edges()

m_teo = p*n*(n-1)/2 # Valor que se espera para la cantidad de enlaces en una red del tipo E-R
    
print('La cantidad m de enlaces difiere del valor esperado en un', np.around(abs(m - m_teo)/m_teo * 100, 4), ' %')


'''
La relación m = p*n*(n-1)/2 es esperada ya que la probabilidad p, para un número grande de nodos, indica la fracción de
nodos que están enlazados respecto del total de pares posibles. Por su parte, el total de pares posibles corresponde al 
combinatorio (n 2) = n*(n-1)/2, el cual indica la cantidad de grupos de dos nodos que se puede formar en un total de n
nodos.
'''

# ii.


'''
Dado que el grado medio <k> de la red es el promedio de los grados de todos los nodos, se obtiene de forma inmediata la 
relación <k> = 1/n * sum(k_i) = 2*m/n. El factor 2 indica que al sumar todos los grados, inevitablemente contamos 2 veces
cada enlace (ya que la red es no dirigida).

'''

k_med = 2*m / n # Grado medio de la red

k_med_teo = p * (n-1) # Valor esperado para el grado medio en una red del tipo E-R

print('El grado medio <k> difiere del valor esperado en un', np.around(abs(k_med - k_med_teo)/k_med_teo * 100, 4), '%')

'''
Vemos que esta relación se cumple como consecuencia inmediata de la utilizada en i.
Si la reemplazamos en el cálculo de <k> nos queda 2* [p*n*(n-1)/2] / n, con lo cual, luego de simplificar, obtenemos la
relación <k> = p * (n-1)
'''

#%%

red_er_lib = nx.erdos_renyi_graph(int(n), p) # Creamos una red E-R utilizando la librería networkx

# i.

m =  red_er_lib.number_of_edges()

m_teo = p*n*(n-1)/2 # Valor que se espera para la cantidad de enlaces en una red del tipo E-R
    
print('La cantidad m de enlaces difiere del valor esperado en un', np.around(abs(m - m_teo)/m_teo * 100, 4), ' %')


k_med = 2*m / n # Grado medio de la red

k_med_teo = p * (n-1) # Valor esperado para el grado medio en una red del tipo E-R

print('El grado medio <k> difiere del valor esperado en un', np.around(abs(k_med - k_med_teo)/k_med_teo * 100, 4), '%')


'''
Vemos que, al igual que con nuestro código, la red E-R generada mediante la librería de networkx cumple las mismas
para m y <k>.
'''

#%%


'''
b)
'''















#%%

################################################################################
#                               PUNTO 3 
################################################################################