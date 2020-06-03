# -*- coding: utf-8 -*-
"""
Created on Fri May 29 08:41:15 2020

@author: sanz1
"""
import networkx as nx
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random
from tqdm import tqdm

#%%

'''
Funciones útiles
'''


def abrir_txt(nombre_archivo):
    archivo=open(nombre_archivo)
    data=[]
    for linea in archivo:
        linea=linea.strip()
        columna=linea.split()
        data.append(columna)    
    return data

def abrir_reguly(nombre_archivo):
    archivo = open(nombre_archivo)
    data = []
    for linea in archivo:
        linea=linea.strip()
        columna=linea.split()
        data.append(columna)  
    data = data[1:]
    enlaces = []
    for linea in data:
      enlaces.append(linea[:2])
    return(enlaces)

def abrir_esenciales(nombre_archivo):
    archivo = open(nombre_archivo)
    data = []
    for linea in archivo:
        linea=linea.strip()
        columna=linea.split()
        data.append(columna)  
    data = data[2:]
    esenciales = []
    for linea in data:
      if len(linea)>2:
        esenciales.append(linea[1])
    return(esenciales)


#%%

#ruta donde se encuentran los archivos descargados:
path = 'C:/Users/Mati/Documents/GitHub/Redes/tc3 data/'

save_path='C:/Users/Mati/Documents/GitHub/Redes/tc3 data/Resultados/'

#lista con los nombres de los archivos:
filename = ['yeast_Y2H','yeast_AP-MS','yeast_LIT','yeast_LIT_Reguly']

grafos = []  #Inicializamos una lista que contendrá los 4 grafos de networkx en cada índice.

for file in filename:
    
    nombre_archivo = path + file + '.txt'
    
    if file!='yeast_LIT_Reguly':
        
        data = abrir_txt(nombre_archivo) #creamos una lista con los enlaces para cada red.
        
    else:
        
        data = abrir_reguly(nombre_archivo)
        
    grafo = nx.Graph() #Inicializamos un grafo en cada paso de la iteración
        
    grafo.add_edges_from(data) #agregamos los nodos y enlaces
    
    grafos.append( grafo ) #guardamos el grafo en un índice de la lista grafos

#%%
'''
Tabla 1 de Zotenko
'''
    
tabla1=pd.DataFrame(index=filename,columns=['Nodos','Enlaces','<k>','<C>'])

c=0

for g in grafos:
    
    grados=[j for i,j in g.degree()]
    
    tabla1['Nodos'][filename[c]]=g.number_of_nodes()
    
    tabla1['Enlaces'][filename[c]]=g.number_of_edges()
    
    tabla1['<k>'][filename[c]]=np.mean(grados)
    
    tabla1['<C>'][filename[c]]=nx.average_clustering(g)
    
    c=c+1
    
tabla1.to_csv(save_path+'tabla_1.csv', encoding='utf-8')

#%%
'''
Tabla 1 de Zotenko
'''


tabla2=pd.DataFrame(index=filename, columns=filename)

#para ver los enlaces en común tenemos en cuenta que en una base de datos puede aparecer
#(i,j) y en la otro (j,i).

for i in np.arange(0,len(grafos)):
    for j in np.arange(i+1,len(grafos)):
        c=0
        for n in grafos[i].edges:
            if n in grafos[j].edges:
                c=c+1
            elif n[::-1] in grafos[j].edges:
                c=c+1
        tabla2[filename[j]][filename[i]]=c/tabla1['Enlaces'][filename[i]]
        tabla2[filename[i]][filename[j]]=c/tabla1['Enlaces'][filename[j]]

tabla2.columns = list(np.arange(len(filename)))
tabla2.reset_index(drop=True,inplace=True)
for i in tabla2.index:
    tabla2[i][i]=filename[i]

tabla2.to_csv(save_path+'tabla_2.csv',encoding='utf-8')

#%%
'''
Figura 1 de Zotenko
'''

#primero importamos las proteínas esenciales

esenciales = abrir_esenciales(path+'Essential_ORFs_paperHe.txt')

#Notas: en un momento divide por 0 eso hay que arregrarlo, tmb habría que iterar
#más de 1 vez por si algunos hubs tienen igual grado y achicar el paso del
#cutoff a 0.0001

df=pd.DataFrame(index=filename,columns=['fraccion','cutoff'])
cutoff=np.arange(0,1,0.01)
h=0
for g in grafos:
    grados=sorted(g.degree, key=lambda x: x[1], reverse=True)
    hubs=[]
    es=[]
    for i in cutoff:
        hubs_count=int(round(i*len(grados)))
        hubs.append(hubs_count)
        grados_aux=grados[0:hubs_count]
        c=0
        for n,m in grados_aux:
            if n in esenciales:
                c=c+1
        es.append(c)
    df['fraccion'][filename[h]]=np.divide(es,hubs)
    df['cutoff'][filename[h]]=cutoff
    h=h+1

plt.figure()
for i in df.index:
    plt.plot(df['cutoff'][i],df['fraccion'][i],label=i)
plt.legend()
plt.ylim((0,1))
plt.show()
plt.close()


#%%

'''
Figura 3 de Zotenko
'''

#grafos = [grafos[0]]


n_romper = 1 # Hasta cuantos nodos queremos quedarnos tras romper la red

# Armamos una lista con las diferentes funciones de centralidad:
centralidades = [nx.eigenvector_centrality, 'grado', 'random', nx.betweenness_centrality, nx.closeness_centrality]
#centralidades = [nx.eigenvector_centrality, nx.betweenness_centrality]
centralidades_str = ['autoval', 'grado', 'random', 'interm', 'cercania']

centralidades_str = ['autoval', 'betweenness']

data_gc = []

nombre_grafo = 0



for grafo in grafos:
    
    print('Trabajando en una nueva red')
    
    tamaños_gc = []
    
    nombre_centr = 0
    
    for centralidad in tqdm(centralidades):
        
        #print('Utilizando nueva centralidad')
        
        gc = grafo.copy() # Creamos una copia para no desarmar la red original
        
        tamaño_gc = []
        
        while len(gc.nodes()) > n_romper:
    
            gc_nodes = max(nx.connected_components(gc), key = len) # Nodos de la componente gigante
            
            nodos_no_gc = set(gc.nodes()) - set(gc_nodes)
            
            gc.remove_nodes_from(nodos_no_gc) # Creamos el subgrafo de la componente gigante
            
            if len(gc.nodes()) != len(gc_nodes):
                
                raise ValueError('Los nodos de la CG no se eliminaron correctamente')
            
            tamaño_gc.append(len(gc.nodes()))
            
            if centralidad == 'grado':
                    
                max_centralidad = sorted(dict(gc.degree()).items(), key=lambda x: x[1], reverse=True)[0][0] # Tomamos el nodo de mayor grado
                
            elif centralidad == 'random':
                
                max_centralidad = random.choice(list(gc.nodes())) # Elegimos cualquier nodo al azar
                
            elif centralidad == nx.eigenvector_centrality:
                
                # Buscamos el nodo con el nodo de mayor centralidad:
                max_centralidad = sorted(centralidad(gc, max_iter = 10000).items(), key=lambda x: x[1], reverse=True)[0][0]
            
            else: 
                
                # Buscamos el nodo con el nodo de mayor centralidad:
                max_centralidad = sorted(centralidad(gc).items(), key=lambda x: x[1], reverse=True)[0][0]
                
            
                    
            gc.remove_node(max_centralidad) # Eliminamos el nodo más central 
        
        np.savetxt(save_path+filename[nombre_grafo]+'_'+centralidades_str[nombre_centr]+'.txt', np.asarray(tamaño_gc))
        
        nombre_centr += 1
        
        tamaños_gc.append(tamaño_gc)
    
    nombre_grafo += 1    
    
    data_gc.append(tamaños_gc)
            
            
          
    
'''
NOTA: el eigenvector value no converge.

    hay 2 modos que siempre dan igual para las dos centralidades diferentes.
    
    Idea:
        separar por cada centralidad. Algo raro pasa al hacer todo junto



'''


#%%

y2h, apms, lit, lit_reg = data_gc
#y2h = data_gc[0]
cont = 0

plt.figure(1)

for i in y2h:
   
    i = np.asarray(i) 
   
    x = np.linspace(1, max(i), len(i))/max(i)
    
    plt.plot(x, i/max(i), label = centralidades_str[cont])
    
    plt.title('yeast_Y2H')
    
    plt.xlabel('fracción de nodos removidos')
    
    plt.ylabel('fracción del tamaño original de la C-G')
    
    plt.legend()

    cont += 1


cont = 0

plt.figure(2)

for i in apms:
    
    i = np.asarray(i)
    
    x = np.linspace(1, max(i), len(i))/max(i)
    
    plt.plot(x, np.asarray(i)/max(i), label = centralidades_str[cont])
    
    plt.title('yeast_AP-MS')
    
    plt.xlabel('fracción de nodos removidos')
    
    plt.ylabel('fracción del tamaño original de la C-G')
    
    plt.legend()

    cont += 1

cont = 0

plt.figure(3)

for i in lit:
    
    i = np.asarray(i)

    x = np.linspace(1, max(i), len(i))/max(i)
    
    plt.plot(x, np.asarray(i)/max(i), label = centralidades_str[cont])
    
    plt.title('yeast_LIT')
    
    plt.xlabel('fracción de nodos removidos')
    
    plt.ylabel('fracción del tamaño original de la C-G')
    
    plt.legend()

    cont += 1

plt.figure(4)

cont = 0

for i in lit_reg:
    
    x = np.linspace(1, max(i), len(i))/max(i)
    
    i = np.asarray(i)
    
    plt.plot(x, np.asarray(i)/max(i), label = centralidades_str[cont])
    
    plt.title('yeast_LIT_Reguly')
    
    plt.xlabel('fracción de nodos removidos')
    
    plt.ylabel('fracción del tamaño original de la C-G')
    
    plt.legend()

    cont += 1


'''
Tabla 3 de Zotenko
'''


#%%

'''
Figura 2b de He
'''








'''
Tabla 5 de Zotenko
'''













