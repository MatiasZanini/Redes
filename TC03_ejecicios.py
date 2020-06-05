# -*- coding: utf-8 -*-
"""
                **Introducción a Redes Complejas en Biología de Sistemas**
                        Trabajo Computacional 3 (Entrega 05/06)


Grupo: Camila Sanz y Matías Zanini.
"""
################################################################################
#                                 PAQUETES 
################################################################################
import networkx as nx
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random
from tqdm import tqdm
from collections import Counter
from sklearn.linear_model import LinearRegression
from scipy.stats import pearsonr


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


#ruta donde se encuentran los archivos descargados:
path = 'D:/Redes 2020/TC03/data/'

save_path='D:/Redes 2020/TC03/'

#nombres de las redes (solo por estética)
names=['Y2H','AP-MS','LIT','LIT Reguly']

#lista con los nombres de los archivos:
filename = ['yeast_Y2H','yeast_AP-MS','yeast_LIT','yeast_LIT_Reguly']

grafos = []  #Inicializamos una lista que contendrá los 4 grafos de networkx en cada índice.

for file in filename:
    
    nombre_archivo = path + file + '.txt'
    
    if file!='yeast_LIT_Reguly':
        
        data = abrir_txt(nombre_archivo) #creamos una lista con los enlaces para cada red.
        
    else:
        
        data=abrir_reguly(nombre_archivo)
        
    grafo = nx.Graph() #Inicializamos un grafo en cada paso de la iteración
        
    grafo.add_edges_from(data) #agregamos los nodos y enlaces
    
    grafos.append( grafo ) #guardamos el grafo en un índice de la lista grafos

#%%
'''
Tabla 2 de Zotenko
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

#Para ver los enlaces en común, tenemos en cuenta que en una base de datos 
#puede aparecer (i,j) y en la otra (j,i).

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

df=pd.DataFrame(index=filename,columns=['fraccion','cutoff'])

cutoff=np.arange(0,1,0.001)

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

h=0

fig, ax = plt.subplots()

for i in df.index:
    
    plt.plot(df['cutoff'][i],df['fraccion'][i],label=names[h])
    
    h=h+1

plt.legend(fontsize=12)

plt.ylim((0,1))

plt.xlim((0,1))

plt.ylabel('Fracción de nodos esenciales',fontsize=15)

plt.xlabel('Fracción de hubs',fontsize=15)

plt.yticks(np.arange(0.2, 1.2, step=0.2))

plt.grid()

ax.tick_params(axis='both', labelsize=15)

plt.savefig(save_path+'Figura 1 Zotenko.png',bbox_inches='tight')

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
    
# centralidades = [nx.closeness_centrality]

# centralidades_str = ['cercania']

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
    
    # Reemplazamos filename[i] con i el numero del archivo que quiero abrir. Cambiarlo pisara los anteriores
    
    y2h.append(np.loadtxt(save_path+filename[3]+'_'+centralidades_str[nombre_centr]+'.txt'))   

for nombre_centr in range(len(centralidades_str)):
        
    apms.append(np.loadtxt(save_path+filename[3]+'_'+centralidades_str[nombre_centr]+'.txt')) 

for nombre_centr in range(len(centralidades_str)):
    
    lit.append(np.loadtxt(save_path+filename[3]+'_'+centralidades_str[nombre_centr]+'.txt')) 

for nombre_centr in range(len(centralidades_str)):
        
    lit_reg.append(np.loadtxt(save_path+filename[3]+'_'+centralidades_str[nombre_centr]+'.txt'))     
    

    
    

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

plt.grid()
#%%
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

plt.grid()

#%%
cont = 0

plt.figure(3)

gc_esencial = ruptura_esencial[2]

plt.grid()

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

#%%
plt.figure(4)

gc_esencial = ruptura_esencial[3]

plt.plot(gc_esencial[0], gc_esencial[1], 'p', markersize=15, label = 'Escenciales')

cont = 0

plt.grid()

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

#%%
al_azar = []

for grafo in grafos:

    grafo_original = grafo.copy()
        
    gc_componentes = sorted(nx.connected_components(grafo_original), key=len, reverse=True)
            
    gc_original = grafo_original.subgraph(gc_componentes[0]).copy() # Esta es la componente gigante original
                    
    N_gc_original = len(gc_original.nodes())

    esenciales_gc = []
    
    for nodo_es in esenciales:
        
        if nodo_es in list(gc_original.nodes()):
            
            esenciales_gc.append(nodo_es)

    esenciales_grados = sorted(gc_original.degree(esenciales_gc), key=lambda x: x[1], reverse=False)

    esenciales_grados_val = [grado for (nodo, grado) in esenciales_grados]
    
    esenciales_grados_nodo = [nodo for (nodo, grado) in esenciales_grados]
    
    bineo = np.arange(0,max(esenciales_grados_val)+1)
    
    distrib, bineo = np.histogram(esenciales_grados_val, bins = bineo)
    
    # distrib contiene la cantidad de nodos esenciales de grado i en su componente i-esima, salvo la ultima. 
    # La ultima se encuentra corrida de forma que distrib[len(distrib)] = distrib[len(distrib)-1]
    
    grados_a_borrar = list( np.where(distrib!=0)[0]) # Lista con los grados que tenemos que eliminar
    
    listas_grados_random = []
    
    nodos_no_esenciales = set(gc_original.nodes()) - set(esenciales_grados_nodo)
    
    grados_no_esenciales = list(gc_original.degree(nodos_no_esenciales)) # Grados de los nodos no esenciales
    
    for i in grados_a_borrar:
        
        if i == max(grados_a_borrar): # Solucionamos el corrimiento del ultimo valor en distrib
            indice_para_distrib = i-1
        
        else:
            
            indice_para_distrib = i
        
        
        lista_grado_i = []
        
        for grado in grados_no_esenciales:
            
            if grado[1] == i:
                
                lista_grado_i.append(grado[0])
        
        if len(lista_grado_i) < distrib[indice_para_distrib]:  # Si tengo menos nodos no esenciales de grado i que esenciales, relleno con esenciales
            
            index = list(np.where(np.asarray(esenciales_grados_val)==i)[0])
            
            for j in index:
                
                while len(lista_grado_i) < distrib[indice_para_distrib]:
                
                # Rellenamos los nodos no esenciales que faltan con nodos esenciales:
                    lista_grado_i.append(esenciales_grados_nodo[j]) 
            
        listas_grados_random.append(lista_grado_i) 
        # Cada elemento de listas_grados_random contiene una lista con los nodos no esenciales de grado igual a los
        # nodos esenciales. En caso de que la lista sea menor a la cantidad que se neecista por la distribucion de grado,
        # se rellena con nodos esenciales, los cuales podran ser borrados "como si fueran no esenciales"
    
    iteraciones = 50 # Cantidad de veces que queremos iterar para luego promediar
    
    tamaños_gc = []
    
    for iteracion in range(iteraciones):
    
        gc = gc_original.copy()
        
        nodos_a_borrar = []
        
        numero_grado = 0
        
        for borrable in listas_grados_random:
            
            if grados_a_borrar[numero_grado] == max(grados_a_borrar):
                
                indice_para_distrib = grados_a_borrar[numero_grado] -1
        
            else:
            
                indice_para_distrib = grados_a_borrar[numero_grado]
            
            cantidad_a_borrar = distrib[indice_para_distrib] # Cuantos nodos del grado del paso actual hay que borrar
            
            a_borrar = random.sample(borrable, cantidad_a_borrar)
            
            nodos_a_borrar = nodos_a_borrar + a_borrar
            
            numero_grado +=1
        
        gc.remove_nodes_from(nodos_a_borrar)
        
        gc_nodes = max(nx.connected_components(gc), key = len) # Nodos de la nueva componente gigante
                            
        nodos_no_gc = set(gc.nodes()) - set(gc_nodes)
                            
        gc.remove_nodes_from(nodos_no_gc) # Creamos el subgrafo de la componente gigante
                           
        if len(gc.nodes()) != len(gc_nodes):
                               
            raise ValueError('Los nodos de la CG no se eliminaron correctamente')
        
        tamaño_gc = len(gc.nodes())/N_gc_original
        
        tamaños_gc.append(tamaño_gc)
    
    tamaño_random = np.mean(tamaños_gc)
    
    desviacion = np.std(tamaños_gc)
    
    al_azar.append([[tamaño_random], [desviacion]])
    
    
for i in range(len(filename)):
    
    tabla_esenciales['Al azar'][filename[i]] = str(( float(al_azar[i][0][0]) )) + ' \pm ' + str( float(al_azar[i][1][0]) )

print(tabla_esenciales)

#%%

#Para pasar la tabla a latex:

    
print(tabla_esenciales.to_latex())
    



#%%

'''
Figura 2b de He
'''
model = LinearRegression()#para hacer el ajuste lineal

parametros=pd.DataFrame(index=filename, columns=['alpha','beta','p'])

h=0

k_max=[8,9,9,9]

for g in grafos:
    
    grados=sorted(g.degree, key=lambda x: x[1], reverse=False)
    
    dict_g=dict(grados)
    
    grados_val=[jj for ii,jj in grados if jj<=k_max[h]]
    
    frame=pd.DataFrame(index=np.arange(0,1),columns=sorted(set(grados_val),reverse=False))
    
    for i in frame.columns:
        
        nodes=[k for k,v in dict_g.items() if v == i]
        
        c=0
        
        for j in nodes:
            
            if j in esenciales:
                
                c=c+1
        
        frame[i][0]=np.log(1-c/len(nodes))
    
    x=np.array(frame.columns)
    
    ajuste=model.fit(x.reshape(-1, 1), np.array(frame.loc[0]))
    
    r_sq = model.score(x.reshape(-1, 1), np.array(frame.loc[0]))
    
    a=ajuste.coef_[0]
    
    b=ajuste.intercept_
    
    y=a*x+b
    
    fig, ax = plt.subplots()
    
    plt.scatter(x,frame.loc[0])
    
    plt.plot(x,y,'r--',label='Fit: a='+str(a)[0:6]+' b='+str(b)[0:6]+', R^2='+str(r_sq)[0:5])
    
    plt.title(names[h],fontsize=20)
    
    plt.ylabel('Ln(1-PE)', fontsize=15)
    
    plt.xlabel('K',fontsize=15)
    
    plt.legend()
    
    plt.grid()

    ax.tick_params(axis='both', labelsize=15)
    
    plt.savefig(save_path+filename[h]+'He 2b.png',bbox_inches='tight')
    
    plt.show()
    
    plt.close()
    
    rcorr, pvalue = pearsonr(frame.loc[0], y)
    
    print(rcorr,pvalue)
    
    alpha=1-np.exp(a)
    
    beta=1-np.exp(b)
    
    print('alpha='+str(alpha)+', beta='+str(beta))
    
    parametros['alpha'][filename[h]]=alpha
    
    parametros['beta'][filename[h]]=beta
    
    parametros['p'][filename[h]]=pvalue
    
    h=h+1

parametros.to_csv(save_path+'alpha_beta.csv',encoding='utf-8')
#%%

'''
Tabla 5 de Zotenko
'''
tabla5=pd.DataFrame(index=filename,columns=['Pares','Pares del mismo tipo','Pares estimados','std'])

std_he=[]

h=0

vecinos=[2,3,3,3]

for g in grafos:
    
    print(h)
    
    nodes=list(g.nodes())
    
    pairs=0
    
    total_pairs=0
    
    for i in np.arange(0,len(nodes)):
        
        for j in np.arange(i+1,len(nodes)):
            
            if (nodes[i],nodes[j]) not in g.edges():
                
                cn=list(nx.common_neighbors(g, nodes[i], nodes[j]))
                
                if len(cn)>=vecinos[h]:
                    total_pairs=total_pairs+1
                
                if len(cn)>=vecinos[h] and nodes[i] in esenciales and nodes[j] in esenciales:
                    pairs=pairs+1
                
                if len(cn)>=vecinos[h] and nodes[i] not in esenciales and nodes[j] not in esenciales:
                    pairs=pairs+1
    
    tabla5['Pares'][filename[h]]=total_pairs
    
    tabla5['Pares del mismo tipo'][filename[h]]=pairs
    
    print(total_pairs,pairs)
    
    #alpha:probabilidad de que una interacción sea esencial y por ende
    #los nodos involucrados sean esenciales.
    #beta: probabilidad de que un nodo sea esencial por otra causa
    
    x=0
    
    he_model=[]
    
    while x<100:
        
        print(x)
        
        total_nodes=[]
        
        edges=list(g.edges())
        
        en_es=parametros['alpha'][filename[h]]*(len(edges))
        
        enlaces_r=random.sample(range(1, len(edges)), int(en_es))
        
        enlaces_r_list=[]
        
        for n in enlaces_r:
        
            aux=edges[n]
            
            enlaces_r_list.append(aux)
            
            for m in aux:
            
                if m not in total_nodes:
                
                    total_nodes.append(m)
        
        H=g.copy()
        
        H.remove_nodes_from(total_nodes)
        
        nodesH=list(H.nodes())
        
        n_es=parametros['beta'][filename[h]]*(len(nodes))
        
        nodos_r=random.sample(range(1, len(nodesH)), int(n_es))
        
        nodos_r_list=[]
        
        for nn in nodos_r:
        
            auxn=nodesH[nn]
            
            nodos_r_list.append(auxn)
            
            if auxn not in total_nodes:
            
                total_nodes.append(auxn)
            
            else:
            
                print('error')
        
        pairs_he=0
        
        for ii in np.arange(0,len(nodes)):
        
            for jj in np.arange(ii+1,len(nodes)):
            
                if (nodes[ii],nodes[jj]) not in g.edges():
                
                    cn=list(nx.common_neighbors(g, nodes[ii], nodes[jj]))
                    
                    if len(cn)>=vecinos[h] and nodes[ii] in total_nodes and nodes[jj] in total_nodes:
                    
                        pairs_he=pairs_he+1
                    
                    if len(cn)>=vecinos[h] and nodes[ii] not in total_nodes and nodes[jj] not in total_nodes:
                    
                        pairs_he=pairs_he+1
       
        he_model.append(pairs_he)
        
        x=x+1
        
    std_he.append(np.std(he_model)/100)
    
    tabla5['Pares estimados'][filename[h]]=round(np.mean(he_model)) 
    
    tabla5['std'][filename[h]]=round(np.std(he_model)) 
    
    h=h+1

print(tabla5)

tabla5.to_csv(save_path+'tabla_5.csv',encoding='utf-8')
    
#%%












