# -*- coding: utf-8 -*-
"""
Created on Sun Jun 28 13:52:49 2020

@author: sanz1
"""

#proyeccion de la red no dirigida sobre los posts.

#%%
#la idea es la siguiente: Hacemos la red general, la proyectamos y dsp levantamos las características generales
#de la red incluida la distribución de grados.

import networkx as nx
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from networkx.algorithms import bipartite
import community
import matplotlib as mpl
import igraph as ig

#%%
def graph_color(grafo,lista_comunidad,titulo,posi,labels_node,peso_norm):#plotea la red con los colores
    low, *_, high = sorted(lista_comunidad)
    norm = mpl.colors.Normalize(vmin=low, vmax=high, clip=True)
    mapper = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.coolwarm)
    
    nx.draw(grafo,node_color=[mapper.to_rgba(i) for i in lista_comunidad],labels=labels_node,
                              width=peso_norm, pos=posi)
    #nx.draw(grafo,node_color=[mapper.to_rgba(i) for i in lista_comunidad],with_labels=True,
                              #width=peso_norm, pos=posi)
    plt.title(titulo)

def graph_color_lab(grafo,lista_comunidad,titulo,posi,labels_node,peso_norm):#plotea la red con los colores
    low, *_, high = sorted(lista_comunidad)
    norm = mpl.colors.Normalize(vmin=low, vmax=high, clip=True)
    mapper = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.coolwarm)
    
    #nx.draw(grafo,node_color=[mapper.to_rgba(i) for i in lista_comunidad],labels=labels_node,
     #                         width=peso_norm, pos=posi)
    nx.draw(grafo,node_color=[mapper.to_rgba(i) for i in lista_comunidad],with_labels=True,
                              width=peso_norm, pos=posi)
    plt.title(titulo)
#%%
path = 'D:/Redes 2020/tp_final_pruebas/Pruebas 2/'

pathF='D:/Redes 2020/tp_final_pruebas/Pruebas 2/red_posts/'

guardar=pd.read_pickle(path+'datos_red.p')

likes=[]
for i in guardar.index:
    likes.append(len(guardar['reacters'][i]))
guardar['likes']=likes

colores_cat=['black','grey','limegreen','magenta','orange']
#la idea es hacer un grafo dirigido con 2 tipos de nodos,personas y post entonces tenemos poster-->post-->reacters

#red original
enlaces=[]
my_dict={}
for i in guardar.index:
    if (guardar['poster'][i],str(i)) not in enlaces:
        enlaces.append((guardar['poster'][i],str(i)))
        my_dict[(guardar['poster'][i],str(i))]=int(guardar['categoria'][i])
    for j in guardar['reacters'][i]:
        if (str(i),j) not in enlaces:
            enlaces.append((str(i),j))
            my_dict[(str(i),j)]=int(guardar['categoria'][i])

G=nx.Graph()#no hace falta hacerla con multienlace (por si un posteador tmb likeo, lo consideramos como uno solo)
G.add_edges_from(enlaces)
nx.set_edge_attributes(G, my_dict, 'categoria')


post_num=[]
for i in np.arange(0,len(guardar)):
    post_num.append(str(i))

nodes_dict={}
color_nodes=[]
for n in G.nodes:
    if n in post_num:
        nodes_dict[n]=0#0 es un post
        ind=post_num.index(n)
        cat=guardar['categoria'][ind]
        color_nodes.append(colores_cat[int(cat)])
    else:
        nodes_dict[n]=1#1 es una persona
        color_nodes.append(colores_cat[0])
        
nx.set_node_attributes(G, nodes_dict, 'post-persona')

#%%
#hasta acá tenemos lo mismo que en la red completa, pero con una red no dirigida.
#entonces ahora vamos a proyectar sobre los posts
bottom_nodes, top_nodes = bipartite.sets(G)
#proy_personas = bipartite.projected_graph(G_nodir, bottom_nodes,multigraph=True)
proy_post=bipartite.projected_graph(G,top_nodes,multigraph=True)

#ahora lo que hacemos es trasnformar esta red multienlace en una red pesada.

A_post = nx.to_pandas_adjacency(proy_post)
w_max=np.max(np.max(A_post))#hacemos que los pesos vayan entre 0 y 1
enlaces2=[]
peso=[]
for i in A_post.index:
    for j in A_post.columns:
        if A_post[j][i]!=0:
            enlaces2.append((i,j,A_post[j][i]/w_max))
            peso.append(A_post[j][i]/w_max)

nuevo_post=nx.Graph()
nuevo_post.add_weighted_edges_from(enlaces2)

color_list=[]#cambiamos y en lugar de ponerle color a los enlaces le ponemos color a los posts
#segun corresponda
my_dict2={}
for n in nuevo_post.nodes():
    extra=int(guardar['categoria'][int(n)])
    color_list.append(colores_cat[extra])
    my_dict2[n]=extra

nx.set_node_attributes(nuevo_post, my_dict2, 'categoria')


#%%
#dibujamos la red

pos = nx.spring_layout(nuevo_post,weight='weight')#yo creo que se ve mucho mejor así
pos2 = nx.kamada_kawai_layout(nuevo_post,weight='weight')#yo creo que se ve mucho mejor así
fig = plt.figure(figsize=[9,9])
nx.draw(nuevo_post,labels=my_dict2,
        width=peso, # Ancho de los enlaces
       node_color = color_list, # Color de los nodos
       edge_color = 'black', # Color de los enlaces
       node_size = 400,pos=pos2)

df=pd.DataFrame(index=['nodos','enlaces','k medio','k max','k min','w medio','w max','w min','densidad'],columns=['value'])
df['value']['nodos']=nuevo_post.number_of_nodes()
df['value']['enlaces']=nuevo_post.number_of_edges()
grados=nuevo_post.degree()
grados_val=[j for i,j in grados]
df['value']['k medio']=np.mean(grados_val)
df['value']['k max']=np.max(grados_val)
df['value']['k min']=np.min(grados_val)
df['value']['densidad']=nx.density(nuevo_post)
df['value']['w medio']=np.mean(peso)
df['value']['w max']=np.max(peso)
df['value']['w min']=np.min(peso)

df.to_excel(pathF+'datos_red_posts.xlsx')
#%%
#calculo la modularidad tomando la partición original
#la formula es la misma pero con los pesos
#no me sale clacularla así que voy a usar la de networkx (hablar con Fede)

#me asombra que no den igual
Q_multi=community.modularity(my_dict2,proy_post)#-0.25256714817011905
Q_peso=community.modularity(my_dict2,nuevo_post)#0.020189548148989732

#qué pasa si ponemos un treshhold



#%%

#aplicamos algoritmos de comunidades

treshold=np.mean(peso)
nuevo_post2=nuevo_post.copy()
peso2=[]
for i in nuevo_post.edges():
    if nuevo_post.edges[i]['weight']<=treshold:
        nuevo_post2.remove_edge(i[0],i[1])
    else:
        peso2.append(nuevo_post.edges[i]['weight'])

v=[]
for i in nuevo_post.nodes():
    if nuevo_post2.degree(i)==0:
        nuevo_post2.remove_node(i)
        v.append((i,guardar['categoria'][int(i)]))
print(nuevo_post.number_of_nodes(),nuevo_post2.number_of_nodes())

#saco 4 de serie y 2 de interno

color_list=[]#cambiamos y en lugar de ponerle color a los enlaces le ponemos color a los posts
#segun corresponda
my_dict2={}
for n in nuevo_post2.nodes():
    extra=int(guardar['categoria'][int(n)])
    color_list.append(colores_cat[extra])
    my_dict2[n]=extra
    
pos = nx.kamada_kawai_layout(nuevo_post2,weight='weight')
fig = plt.figure(figsize=[9,9])
nx.draw(nuevo_post2,labels=my_dict2,
        width=peso2, # Ancho de los enlaces
       node_color = color_list, # Color de los nodos
       edge_color = 'black', # Color de los enlaces
       node_size = 400,pos=pos)


Q_tresh=community.modularity(my_dict2,nuevo_post2,weight='weight')#0.04

#%%

#después de pensarlo un rato, me parece que kmeans sería un buen algoritmo para usar

#para mi se ve muy bien con louvain y tenemos que probar qué onda con kmeans, o con eb-infomap
#imponiendo el corte del dendograma.
louvain=community.best_partition(nuevo_post2,weight='weight')
lista_com_lou=list(louvain.values())
print('Louvain') 
print('Cantidad de comunidades óptimas: '+str(np.max(lista_com_lou)+1))
fig = plt.figure(figsize=[9,9])
graph_color(nuevo_post2,lista_com_lou,'Louvain',pos,my_dict2,peso2)
Q_l=community.modularity(louvain,nuevo_post2,weight='weight')#0.071
fig = plt.figure(figsize=[9,9])
graph_color_lab(nuevo_post2,lista_com_lou,'Louvain',pos,my_dict2,peso2)



post_ig = ig.Graph.TupleList(edges=nuevo_post2.edges(),directed=False)
edges=ig.EdgeSeq(post_ig)
vseq = post_ig.vs #lista de nodos (vseq['name'])
pesos_ig=[]

for i in edges:
    aux=(vseq['name'][i.tuple[0]],vseq['name'][i.tuple[1]])
    pesos_ig.append(nuevo_post2.edges[aux]['weight'])

post_ig.es['weight'] = pesos_ig

def nodos_nx_ig(nodos_nx,nodos_ig,comunidad):#paso de los nodos ordenados segun ig a los nodos ordenados segun nx
    lista=[None]*len(nodos_nx)
    for i in np.arange(0,len(comunidad)):
        for j in comunidad[i]:
            node_name=nodos_ig[j]['name']
            ind=nodos_nx.index(node_name)
            lista[ind]=i
    return lista

nodos=list(nuevo_post2.nodes())

fg_dendograma=post_ig.community_fastgreedy(weights=post_ig.es['weight'])
print('Fast Greedy') 
print('Cantidad de comunidades óptimas: '+str(fg_dendograma.optimal_count))
comunidades_fg=fg_dendograma.as_clustering()#comunidades
lista_com_fg=nodos_nx_ig(nodos,vseq,comunidades_fg)
fig = plt.figure(figsize=[9,9])
graph_color(nuevo_post2,lista_com_fg,'Fast Greedy',pos,my_dict2,peso2)

eb_dendograma = post_ig.community_edge_betweenness(clusters=None, directed=False, weights=post_ig.es['weight'])
print('Edge betweenness')
print('Cantidad de comunidades óptimas: '+str(eb_dendograma.optimal_count))
comunidades_eb=eb_dendograma.as_clustering()
lista_com_eb=nodos_nx_ig(nodos,vseq,comunidades_eb)
fig = plt.figure(figsize=[9,9])
graph_color(nuevo_post2,lista_com_eb,'EB',pos,my_dict2,peso2)

infomap=post_ig.community_infomap(edge_weights=post_ig.es['weight'], trials=1000)
print('Infomap')
print('Cantidad de comunidades óptimas: '+str(len(infomap)))


post_ig = ig.Graph.TupleList(nuevo_post.edges(), directed=False)
vseq = post_ig.vs #lista de nodos (vseq['name'])

#Louvain
louvain=community_louvain.best_partition(nuevo_post)
lista_com_lou=list(louvain.values())
print('Louvain') 
print('Cantidad de comunidades óptimas: '+str(np.max(lista_com_lou)+1))

#Fast Greedy
fg_dendograma=G_igraph.community_fastgreedy(weights=None)
print('Fast Greedy') 
print('Cantidad de comunidades óptimas: '+str(fg_dendograma.optimal_count))
comunidades_fg=fg_dendograma.as_clustering()#comunidades

lista_com_fg=nodos_nx_ig(nodos,vseq,comunidades_fg)

#Edge beetweeness
eb_dendograma = G_igraph.community_edge_betweenness(clusters=None, directed=False, weights=None)

print('Edge betweenness')

print('Cantidad de comunidades óptimas: '+str(eb_dendograma.optimal_count))

graficar_dendograma(eb_dendograma)

comunidades_eb=eb_dendograma.as_clustering()

lista_com_eb=nodos_nx_ig(nodos,vseq,comunidades_eb)

#Infomap
infomap=G_igraph.community_infomap()

print('Infomap')

print('Cantidad de comunidades óptimas: '+str(len(infomap)))

lista_com_info=nodos_nx_ig(nodos,vseq,infomap)

#figuras
fig = plt.figure(figsize=[15,15])
plt.subplot(2, 2, 1)
graph_color(G,lista_com_lou,'Louvain',posiciones)
plt.subplot(2, 2, 2)
graph_color(G,lista_com_fg,'Fast Greedy',posiciones)
plt.subplot(2, 2, 3)
graph_color(G,lista_com_eb,'Edge betweenness',posiciones)
plt.subplot(2, 2, 4)
graph_color(G,lista_com_info,'Infomap',posiciones)
plt.savefig(path+'particiones.png')
plt.show()
plt.close()