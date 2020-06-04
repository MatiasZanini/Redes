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

# Armamos una lista con las diferentes funciones de centralidad:
centralidades = [nx.degree_centrality, nx.eigenvector_centrality, 'random', nx.betweenness_centrality, nx.closeness_centrality]

centralidades_str = ['grado', 'autoval', 'random', 'interm', 'cercania']

#Prueba con solo grado:
    
centralidades = [nx.degree_centrality, 'random']

centralidades_str = ['grado', 'random']

data_completa = []

nombre_grafo = 0 # Contador para nombrar cada red

for grafo in grafos:
    
    print('Trabajando en la red', filename[nombre_grafo])
    
    centralidades_por_grafo = []
    
    nombre_centr = 0 # Contador para nombrar cada analisis de centralidad
    
    for centralidad in tqdm(centralidades):
        
        #print('Utilizando nueva centralidad')
        
        grafo_original = grafo.copy() # Creamos una copia para no desarmar la red original
        
        tamaño_gc = [1] # Paso cero. La fraccion tamaño_cg/tamaño_cg_original es 1
        
        nodos_remov = [0] # Paso cero. La fraccion nodos_removidos/nodos_totales_de_la_cg es 0
        
        gc_componentes = sorted(nx.connected_components(grafo_original), key=len, reverse=True)
        
        gc = grafo_original.subgraph(gc_componentes[0]).copy() # Esta es la componente gigante original
                
        N_gc_original = len(gc.nodes())
        
        paso = 1
        
        while paso < N_gc_original/2:
            
            if paso %10 ==0:
                  
                  print('Removidos', paso, 'nodos de', N_gc_original/2) # Mostramos el progreso de la iteración
            
            if len(gc.nodes()) != 0:
                
                if centralidad == 'random':
                    
                    max_centralidad = random.choice(list(gc.nodes())) # Elegimos cualquier nodo al azar
                    
                elif centralidad == nx.eigenvector_centrality:
                    
                    # Buscamos el nodo con mayor centralidad:
                    max_centralidad = sorted(centralidad(gc, max_iter = 10000).items(), key=lambda x: x[1], reverse=True)[0][0]
                
                else: 
                    
                    # Buscamos el nodo con mayor centralidad:
                    max_centralidad = sorted(centralidad(gc).items(), key=lambda x: x[1], reverse=True)[0][0]
                    
                gc.remove_node(max_centralidad) # Eliminamos el nodo más central 
                
                if len(gc.nodes()) !=0:
                
                    conectado = nx.is_connected(gc) # Chequeamos si rompimos la red. 
                    
                    if not conectado: # Si la red sigue conectada no es necesario crear un subgrafo.
                        
                        # Recalculamos la componente gigante:
                        gc_nodes = max(nx.connected_components(gc), key = len) # Nodos de la componente gigante
                        
                        nodos_no_gc = set(gc.nodes()) - set(gc_nodes)
                        
                        gc.remove_nodes_from(nodos_no_gc) # Creamos el subgrafo de la componente gigante
                        
                        if len(gc.nodes()) != len(gc_nodes):
                            
                            raise ValueError('Los nodos de la CG no se eliminaron correctamente')
            
            tamaño_gc.append(len(gc.nodes()) / N_gc_original )
            
            nodos_remov.append(paso / N_gc_original)
            
            paso += 1
           
        
        x_data = np.asarray(nodos_remov)
        
        y_data = np.asarray(tamaño_gc)
        
        data_centralidad = np.asarray([x_data, y_data])
        
        np.savetxt(save_path+filename[nombre_grafo]+'_'+centralidades_str[nombre_centr]+'.txt', data_centralidad)
        
        nombre_centr += 1
        
        centralidades_por_grafo.append(data_centralidad)
    
    nombre_grafo += 1    
    
    data_completa.append(centralidades_por_grafo)
    
    
#%%

# Eliminamos los escenciales:
    
esenciales = abrir_esenciales(path+'Essential_ORFs_paperHe.txt')

nombre_grafo = 0 # Contador para nombrar cada red

ruptura_esencial = []

for grafo in grafos:      
    
    grafo_original = grafo.copy()
    
    gc_componentes = sorted(nx.connected_components(grafo_original), key=len, reverse=True)
        
    gc = grafo_original.subgraph(gc_componentes[0]).copy() # Esta es la componente gigante original
                
    N_gc_original = len(gc.nodes())
    
    nodos_borrar = []
    
    for nodo_es in esenciales:
        
        if nodo_es in list(gc.nodes()):
            
            nodos_borrar.append(nodo_es)
    
    gc.remove_nodes_from(nodos_borrar) # Borramos los nodos esenciales de la componente gigante
     
    gc_nodes = max(nx.connected_components(gc), key = len) # Nodos de la nueva componente gigante
                        
    nodos_no_gc = set(gc.nodes()) - set(gc_nodes)
                        
    gc.remove_nodes_from(nodos_no_gc) # Creamos el subgrafo de la componente gigante
                       
    if len(gc.nodes()) != len(gc_nodes):
                           
        raise ValueError('Los nodos de la CG no se eliminaron correctamente')
        
    y_data = len(gc.nodes())/N_gc_original
    
    x_data = len(nodos_borrar)/N_gc_original
    
    data_esencial = np.array([x_data, y_data])
    
    np.savetxt(save_path+filename[nombre_grafo]+'_escenciales.txt', data_esencial)
   
    ruptura_esencial.append(data_esencial)
    
    nombre_grafo +=1
    
#%%

# Opcional: cargamos desde los archivos guardados

y2h = []

apms = []

lit = []

lit_reg = []
  
for nombre_centr in range(len(centralidades_str)):
        
    y2h.append(np.loadtxt(save_path+filename[0]+'_'+centralidades_str[nombre_centr]+'.txt'))   

for nombre_centr in range(len(centralidades_str)):
        
    apms.append(np.loadtxt(save_path+filename[0]+'_'+centralidades_str[nombre_centr]+'.txt')) 

for nombre_centr in range(len(centralidades_str)):
    
    lit.append(np.loadtxt(save_path+filename[0]+'_'+centralidades_str[nombre_centr]+'.txt')) 

for nombre_centr in range(len(centralidades_str)):
        
    lit_reg.append(np.loadtxt(save_path+filename[0]+'_'+centralidades_str[nombre_centr]+'.txt'))     
    
    
    
    

#%%

y2h, apms, lit, lit_reg = data_completa

#%%

cont = 0

plt.figure(1)

gc_esencial = ruptura_esencial[0]

plt.plot(gc_esencial[0], gc_esencial[1], 'p', markersize=15, label = 'Escenciales')

for data in y2h:
   
    x = data[0]
    
    y = data[1]
    
    plt.plot(x, y, label = centralidades_str[cont])
    
    plt.title('yeast_Y2H')
    
    plt.xlabel('fracción de nodos removidos')
    
    plt.ylabel('fracción del tamaño original de la C-G')
    
    plt.legend()

    cont += 1


cont = 0

plt.figure(2)

gc_esencial = ruptura_esencial[1]

plt.plot(gc_esencial[0], gc_esencial[1], 'p', markersize=15, label = 'Escenciales')

for data in apms:
    
    x = data[0]
    
    y = data[1]
    
    plt.plot(x, y, label = centralidades_str[cont])
    
    plt.title('yeast_AP-MS')
    
    plt.xlabel('fracción de nodos removidos')
    
    plt.ylabel('fracción del tamaño original de la C-G')
    
    plt.legend()

    cont += 1

cont = 0

plt.figure(3)

gc_esencial = ruptura_esencial[2]

plt.plot(gc_esencial[0], gc_esencial[1], 'p', markersize=15, label = 'Escenciales')

for data in lit:
    
    x = data[0]
    
    y = data[1]
    
    plt.plot(x, y, label = centralidades_str[cont])
    
    plt.title('yeast_LIT')
    
    plt.xlabel('fracción de nodos removidos')
    
    plt.ylabel('fracción del tamaño original de la C-G')
    
    plt.legend()

    cont += 1

plt.figure(4)

gc_esencial = ruptura_esencial[3]

plt.plot(gc_esencial[0], gc_esencial[1], 'p', markersize=15, label = 'Escenciales')

cont = 0

for data in lit_reg:
    
    x = data[0]
    
    y = data[1]
    
    plt.plot(x, y, label = centralidades_str[cont])
    
    plt.title('yeast_LIT_Reguly')
    
    plt.xlabel('fracción de nodos removidos')
    
    plt.ylabel('fracción del tamaño original de la C-G')
    
    plt.legend()

    cont += 1

#%%

'''
Tabla 3 de Zotenko
'''

tabla_esenciales = pd.DataFrame(index = filename, columns=['Nodos esenciales','Al azar'])

for i in range(len(filename)):
    
    tabla_esenciales['Nodos esenciales'][filename[i]] = ruptura_esencial[i][1]



for grafo in grafos:

    grafo_original = grafo.copy()
        
    gc_componentes = sorted(nx.connected_components(grafo_original), key=len, reverse=True)
            
    gc_original = grafo_original.subgraph(gc_componentes[0]).copy() # Esta es la componente gigante original
                    
    N_gc_original = len(gc_original.nodes())

    nodos_borrar = []
    
    for nodo_es in esenciales:
        
        if nodo_es in list(gc.nodes()):
            
            nodos_borrar.append(nodo_es)






#%%

'''
Figura 3 de Zotenko -  Cami
'''

df_grado=pd.DataFrame(index=filename,columns=['x','y'])

h=0

#por grado
for g in grafos:

    #separamos la componente gigante
    copia = g.copy()

    Gcc = sorted(nx.connected_components(copia), key=len, reverse=True)

    original=len(Gcc[0])

    sub_g=copia.subgraph(Gcc[0]).copy()

    #por grado
    y=[1]

    x=[0]

    n=1

    while n<original/2:

        medida = nx.degree_centrality(copia)

        #en medida orden tenemos solo los nombres de la proteína
        medida_orden=[key for key, value in sorted(medida.items(), key=lambda item: item[1], reverse=True)]

        #saco el nodo de mayor centralidad
        copia.remove_node(medida_orden[0])
        
        if medida_orden[0] in list(sub_g.nodes()):
            
            sub_g.remove_node(medida_orden[0])
            
            #recalculo la componente gigente
            Gcc_2 = sorted(nx.connected_components(sub_g), key=len, reverse=True)
            
            #computo de nuevo sub_g
            sub_g=sub_g.subgraph(Gcc_2[0]).copy()
            
            #guardo los resultados
            y.append(len(Gcc_2[0])/original)
        
        else:
        
            y.append(y[n-1])
        
        x.append(n/original)
        
        n=n+1
    
    df_grado['x'][filename[h]]=x
    
    df_grado['y'][filename[h]]=y
    
    h=h+1

df_grado.to_pickle(save_path+'centralidades/centralidad_grado.p')


#%%

'''
Figura 2b de He
'''





#%%
'''
Tabla 5 de Zotenko
'''













