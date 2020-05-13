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
import plfit
from collections import Counter
# Evitar acumulación de mensajes de warning en el display
import warnings  
warnings.filterwarnings("ignore")


#%%

################################################################################
#                               PUNTO 1 
################################################################################
'''
Inciso a)
'''

path = 'D:/Redes 2020/TC02/TC02_data/' #Colocar la ruta donde están guardados los archivos

filename=['as-22july06_edgelist','netscience_edgelist','power_enlaces','yeast_Y2H']

grafos=[]

for file in filename:
    
    with open(path+file+'.txt') as f:#abrimos el archivo
        
        data=f.readlines()
    
    for i in np.arange(len(data)):#transformamos en una lista de tuples
        
        if file==filename[1]:
            
            data[i]=data[i].strip().split(';')#el segundo archivo tiene separación con ";"
        
        else:
            
            data[i]=data[i].strip().split()#el resto tiene separación con espacio 
        
        data[i]=tuple(data[i][0:2])#ignoramos los pesos

    G=nx.Graph()#inicializamos el grafo
    
    G.add_edges_from(data)#agregamos los nodos y enlaces
    
    grafos.append(G)#Guardamos los 4 grafos en una lista

#Distribución de grado P(k) (hacemos un histograma y un scatter plot)

for i in np.arange(len(grafos)):
    
    #histograma:
    
    grados=grafos[i].degree
    
    x_degree=[j for k,j in grados]
    
    x_log=np.log10(x_degree)
    
    logbins = np.logspace(np.min(x_log),np.max(x_log),13)
    
    plt.figure()
    
    hist,bins,bars=plt.hist(x_degree,bins=logbins,density=True,facecolor='blue', alpha=0.5, ec='black')
    
    plt.xscale('log')
    
    plt.yscale('log')
    
    plt.title('Histograma de grado '+filename[i])
    
    plt.xlabel('k (log scale)')
    
    plt.ylabel('P(k) (log scale)')
    
    plt.show()
    
    plt.close()
    
    #scatter plot
    
    N=grafos[i].number_of_nodes()
    
    count=Counter(x_degree)
    
    count_orden=sorted(count.items())
    
    k=[]
    
    p_k=[]
    
    for n,m in count_orden:
        
        k.append(n)
        
        p_k.append(m/N)

    plt.figure()
    
    plt.scatter(k,p_k)
    
    plt.xscale('log')
    
    plt.yscale('log')
    
    plt.ylim((np.min(p_k)-10**-(abs(int(math.log10(np.min(p_k))))+1),np.max(p_k)+1))
    
    plt.title('Distribución de probabilidad '+filename[i])
    
    plt.xlabel('k (log scale)')
    
    plt.ylabel('P(k) (log scale)')
    
    plt.show()
    
    plt.close()
    
#%%
'''
Inciso b)
Cualitativamente, se puede ver que la distribución de grados de la red de Internet (as-22july06) y 
de la red de proteínas (yeast_Y2H) corresponeden a una power law. Para verificar estas observaciones, podemos
hacer los ajustes correspondientes (inciso d)
'''

'''
Inciso c)
Todas las redes exiben efectos de borde, por un lado, ninguna red presenta nodos con k=0. Por otro lado, 
si cualitativamente trazamos una recta lineal en los histogramas partiendo del mínimo grado en donde
podemos comenzar a trazar la recta estimamos que deberíamos ver mayor cantidad de hubs en la red.
Cualitativamente, la que parece indicar un efecto de borde mayor es la red de Internet, seguida por la
de proteínas.
Verificamos esta última estimación en el final del script donde comparamos el K máximo que presenta cada
red con el K máximo estimado para redes que siguen una power law.
'''
#%%
'''
Inciso d)
'''
#Guardamos kminimo y gamma de cada red para ver cuantitavente las estimaciones del ejercicio c).

kminimo=[]

gammas=[]

for i in np.arange(len(grafos)):
    
    grados=grafos[i].degree
    
    x_degree=[j for k,j in grados]
    
    fit=plfit.plfit(x_degree)
    
    plt.figure()
    
    fit.plotpdf()
    
    plt.title('Ajuste '+filename[i])
    
    plt.xlabel('k (log scale)')
    
    plt.ylabel('P(k) (log scale)')
    
    plt.show()
    
    plt.close()
    
    xmin = fit._xmin
    
    alpha = fit._alpha

    print('Red '+filename[i]+': K_min = '+str(xmin)+'; Gamma = '+str(alpha))
    
    gammas.append(alpha)
    
    kminimo.append(xmin)

#%%
'''
Las redes que siguen una power law son invariantes de escala, es decir, para todas las escalas observamos
el mismo comportamiento, esto también implica que no existe una escala característica en el sistema.
Por lo tanto, carece de sentido calcular el grado medio de los nodos de la red ya que las fluctuaciones
de esta variable son de órden de magnitud similar o mayor a la variable misma. 
Para verlo podemos hacer:
'''
for i in [0,3]:#tomamos las 2 redes que siguen una power law
    
    print('Red '+filename[i])
    
    aux=[m for n,m in grafos[i].degree() if m>=kminimo[i]]
    
    promK=np.mean(aux)
    
    promK2=np.mean(np.array(aux)**2)
    
    sigma=(promK2-promK**2)**(1/2)
    
    print('Sigma_k: '+str(sigma))
    
    print('k = '+str(promK)+'+-'+str(sigma))

#%%
#Efectos de bordes en el grado máximo de la red:

for i in np.arange(len(grafos)):
    
    print('Red '+filename[i])
    
    aux=[m for n,m in grafos[i].degree()]
    
    print('K máximo real: '+str(np.max(aux)))
    
    estimado=kminimo[i]*(grafos[i].number_of_nodes()**(1/(gammas[i]-1)))
    
    print('K máximo estimado: '+str(estimado))
   
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