# -*- coding: utf-8 -*-
"""
                **Introducción a Redes Complejas en Biología de Sistemas**
                        Trabajo Computacional 1 (Entrega 05/05)


Grupo: Camila Sanz, Matías Zanini, Debora Copa
"""
################################################################################
#                                 PAQUETES 
################################################################################

import scipy.io as spio   #cargar matlab
import pandas as pd       #DB
import numpy as np        #matemática, simil maltab 
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import igraph as ig
# Evitar acumulación de mensajes de warning en el display
import warnings  
warnings.filterwarnings("ignore")


#%%

################################################################################
#                               PUNTO 1 
################################################################################

def abrir_txt(lista_de_enlaces_txt):

    '''
   Crea una lista con los datos de las redes guardados en .txt
    '''
    
    archivo = open(lista_de_enlaces_txt)
    
    data=[]
    
    for linea in archivo:
    
        linea = linea.strip()
        
        columna = linea.split()
        
        data.append(columna)    
    
    return data


'''
Inciso a)
'''


#ruta donde se encuentran los archivos descargados:
path = 'C:/Users/Mati/Documents/GitHub/Redes/tc1 data/'

#lista con los nombres de los archivos:
filename = ['yeast_AP-MS','yeast_LIT','yeast_Y2H']

grafos = []  #Inicializamos una lista que contendrá los 3 grafos de networkx en cada índice.

for file in filename:
    
    nombre_archivo = path + file + '.txt'
    
    data = abrir_txt(nombre_archivo) #creamos una lista con los enlaces para cada red.
    
    grafo = nx.Graph() #Inicializamos un grafo en cada paso de la iteración
    
    grafo.add_edges_from(data) #agregamos los nodos y enlaces
    
    grafos.append( grafo ) #guardamos el grafo en un índice de la lista grafos



#Ahora graficamos las 3 redes y las comparamos:
   
fig, axis = plt.subplots(1, len(grafos))

contador = 0

for G in grafos:
    
    contador +=1
    
    
    #plt.figure()
   # axis[contador-1].plot(1, len(grafos), contador)
    axis[contador-1].set_title(filename[contador-1], fontsize = 16)
    
    nx.draw(G, 
        width = 2, # Ancho de los enlaces
        node_color = 'blue', # Color de los nodos
        edge_color = 'violet', # Color de los enlaces
        node_size = 10, # Tamaño de los nodos
        with_labels = False, # True o False, Pone o no el nombre de los nodos
        ax = axis[contador-1]
        )

plt.show()    
#Descomentar lo siguiente si se quiere guardar la figura en formato png:
#plt.savefig(path+'1a_comparacion_grafica.png')
    


#%%

'''
Inciso b)
'''

#Ordenamos la tabla poniendo en las columnas los datos que solicita el ejercicio:
df = pd.DataFrame(index = filename, columns = ['N','L','<k>','k_max','k_min','Densidad','<C_i>','C','Diámetro'])

contador = 0

for file in filename:
      
    N = grafos[contador].number_of_nodes() #Numero de nodos
    
    L = grafos[contador].number_of_edges() #Numero de enlaces   
    
    df['N'][file] = N
    
    df['L'][file] = L
    
    df['<k>'][file] = 2*L / N
    
    grados = [m for n,m in grafos[contador].degree()]
    
    df['k_max'][file] = np.max(grados)
    
    df['k_min'][file] = np.min(grados)
    
    df['Densidad'][file] = 2*L  / ( N*(N-1)) #Otra posibilidad es: nx.density(grafos[contador])
    
    df['<C_i>'][file] = nx.average_clustering(grafos[contador])
    
    df['C'][file] = nx.transitivity(grafos[contador])
    
    if nx.is_connected(grafos[contador]):
        df['Diámetro'][file] = nx.diameter(grafos[contador])
    else:
        df['Diámetro'][file] = 'infinito'
    
    contador += 1
    
print(df)

#Descomentar lo siguiente si se desea guardar el archivo con la tabla en formato .csv :
#df.to_csv(path + '1b.csv')

#%%
'''
En cuanto a si la red es dirigida o no, podemos basarnos en la matriz de adyacencia. Si la misma es simétrica, resulta
natural considerar la red como no dirigida, dado que no se tiene en cuenta el sentido de sus enlaces.
Por lo tanto, chequeamos que las 3 redes son simétricas a partir de su matriz de adyacencia:
'''

matrices = []

for G in grafos:
    
    matrices.append(nx.to_numpy_matrix(G))

for i in matrices:
    
    A = i - np.transpose(i) #Si la matriz es igual a su transpuesta, entonces A será nula.

    if np.all(A==0):
        
        print('La matriz es simétrica. Por lo tanto, consideramos la red como no dirigida')
        
    else:
        
        print('La matriz no es simétrica. La red es dirigida')




#%%

'''
Inciso c)

La estructura de las redes es similar en todos los casos. Se encuentra una componente gigante de proteínas de
levadura y varias componentes más pequeñas. 

En terminos de nodos, la red Y2H tiene aproximadamente 500 más que las otras.
Se puede observar que en la red AP-MS, donde los enlaces se adjudican si dos proteínas interactúan entre sí, 
la cantidad de enlaces y el grado medio <k> es superior a los de los demás. Lo mismo ocurre con su grado máximo. 
Por otro lado, la densidad resulta ser similar entre las 3 redes, siendo la AP-MS un poco superior.
Esto era esperable tras observar que posee el mayor grado medio y la mayor cantidad de enlaces. 

El coeficiente de clustering mustra que hay mayor cantidad de nodos que cumplen la clausura transitiva en la AP-MS.
Esto quiere decir que, si la proteína X interactúa con la proteína Y y la Z, es probable que las proteínas Y y Z
interactúen entre sí. Para la red LIT, vemos que el coeficiente de clustering es menor (aproximadamente 1/3 o 1/4 de 
los nodos cumplen clausura transitiva). Una posibilidad es que esto ocurra debido a que una proteína puede
cumplir más de una función. Luego, aquellas proteínas cuyos vecinos cumplen diferentes funciones, no se enlazan
entre sí.
Por último el menor coeficiente de clustering se halla en la red LIT. Esto quiere decir que, si
las proteínas X-Y aparecen en un paper y X-Z aparece en otro, no necesariamente la probabilidad de que 
Y-Z aparezcan en otro paper es grande. Esto puede deberse principalmente a que los papers pueden estar hablando
de diferentes temáticas. Algo similar se puede deducir del coeficiente de clustering global.

En cuanto al diámetro de las redes, dado que las mismas no son conexas, se tomó por convención que su diametro
es infinito.
'''


#%%

'''
Inciso d)
'''

ne = pd.DataFrame(index = filename, columns = ['nodos','enlaces']) #Preparamos un dataframe para los nodos y enlaces

contador = 0

for file in filename:

    ne.loc[file] = [grafos[contador].nodes(), grafos[contador].edges()]
    
    contador += 1

plt.figure()

nodos_venn = venn3([set(ne['nodos'][filename[0]]), set(
    ne['nodos'][filename[1]]),
    set(ne['nodos'][filename[2]])],
    set_labels = ('yeast_AP-MS', 'yeast_LIT', 'yeast_Y2H'))

plt.title('Diagrama de Venn (nodos)')

#Descomentar lo siguiente si se desea guardar la figura como .png :
#plt.savefig(path + '1d_nodos.png')

plt.figure()

nodos_venn = venn3([set(ne['enlaces'][filename[0]]), set(
    ne['enlaces'][filename[1]]),
    set(ne['enlaces'][filename[2]])],
    set_labels = ('yeast_AP-MS', 'yeast_LIT', 'yeast_Y2H'))

plt.title('Diagrama de Venn (enlaces)')

#Descomentar lo siguiente si se desea guardar la figura como .png :
#plt.savefig(path+ '1d_enlaces.png')



#%%

################################################################################
#                               PUNTO 2 
################################################################################



"""#Ejercicio 2
Considere la red social de 62 delfines de Nueva Zelanda (dolphins.txt).

a. Examine diferentes opciones de layout para este grafo e identifique la que le resulte más
informativa. Justifique su elección detallando las características estructurales de la red
que su elección pone en evidencia. Incluya en la representación gráfica de la red información sobre el sexo de los delfines.

b. Se trata una red donde prevalece la homofilia en la variable género? Para responder

i. Considere la distribución nula para la fracción de enlaces que vinculan géneros diferentes, generada a partir de al menos 1000 asignaciones aleatorias de género.

ii. A partir de lo obtenido proponga una estimación para el valor y el error de dicha
cantidad cuando no existe vínculo entre topolgía de la red medio y asignación de
género. Compare su estimación con el valor medio esperado.

iii. Estime la significancia estadística (p-valor) del valor observado en el caso de la
red real.

c. (*) Identifique alguna metodología basada en observables topológicos para eliminar
nodos secuencialmente de la red de manera de dividirla en dos componentes de tamaños
comparables en el menor número de pasos. Explique y muestre los resultados obtenidos.
Intente cuantificar su estrategia comparándola con lo que se obtendría al eliminar nodos
de manera aleatoria.


---------------------------------------------------------

En este ejercicio se busca estudiar una población de delfines y tratar de entender cómo son los vínculos entre dicha población y si dichos vínculos tienen que ver con el sexo de los delfines. Basicamente, queremos entender si la población de delfines es homofílica
##Inciso (a)
La idea de este inciso es explorar las distintas posibilidades de layout que nos otorga networkx en función de determinar cuál es la mejor para dar cuenta de la estructura subyacente a la red.
Existen múltilpes layouts, la idea es que probemos tres o cuatro. Algunos de ellos son: 'spring', 'random', 'circle'.
La nomenclatura para aplicar cada uno de ellos sería:




```
# Habiendo importado las librerías networkx y matplotlib y habiendo generado los grafos basta con:
nx.draw(Red_delfines, layout = 'layout_1')
plt.show()
nx.draw(Red_delfines, layout = 'layout_2')
plt.show()
nx.draw(Red_delfines, layout = 'layout_3')
plt.show()
```

##Inciso (b)
La idea de este inciso es estudiar si existe o no homofilia en la red de delfines. La idea es estudiar la fracción de enlaces, sobre el total, que vincule delfines del mismo sexo. Una posibilidad es contar por separado aquellos enlaces que vinculan macho con macho, con los de hembra con hembra. Lo importante es comprender que resulta necesario saber si este valor es grande o chico, saber con qué compararlo. Para eso es necesario repetir el cálculo anterior (averiguar fracción de enlaces entre mismo sexo) sobre redes aleatorias. Pero enteniendo aleatoriedad en el sentido de romper algunas de las relaciones existentes en la red real pero no todas. Por ejemplo, una posibilidad es asignar los géneros aleatoriamente entre los delfines, utilizando la distribución real de géneros. Otra posibilidad es recablear la red, manteniendo la distirbución de grado intacta. Para esta última alternativa existen funciones en networkx.

##Inciso (c)*
Este inciso es opcional, pero no por ello menos importante (a no preocuparse, que en otros TPs se pide lo mismo). La idea es encontrar una estrategia para romper la red en la menor cantidad de pasos posibles. Es decir, ir eliminando nodos, o enlaces, de manera iterativa, estudiando el tamaño de la componente gigante paso a paso. Las estrategias en este tipo de trabajos se basan en el concepto de centralidad de los nodos, o enlaces, en la red.

# Ejercicio 3
3) Considere la red as-22july06.gml creada por Mark Newman que contiene la estructura de los
sistemas autónomos de internet relevada a mediados de 2006.

a. Encuentre gráficamente la distribución de grado Pk como función de k explorando diferentes alternativas: un bineado lineal o logarítmico, utilizando escalas logarítmicas o lineales en uno o ambos ejes. Discuta que alternativa permite apreciar mejor el carácter libre de escala de dicha distribución.

b. Utilizando funcionalidad de la librería igraph, estime el exponente de dicha distribución

*Comentario Debora:*
*graficar con plt.bar pero que se va a ver mal, mejor con plt.scatter*

"""

import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import random
#from scipy.stats import norm

G = nx.read_gml('D:/Redes 2020/TC01_data/dolphins.gml')
with open('D:/Redes 2020/TC01_data/dolphinsGender.txt') as f:
    gender=f.readlines()
lista_N=[]
lista_G=[]
for i in np.arange(len(gender)):
    gender[i]=gender[i].replace('\n','').split('\t')
    lista_N.append(gender[i][0])
    lista_G.append(gender[i][1])
gender=dict(gender)
plt.figure('Original')
nx.draw(G,labels=gender)
plt.title('Original')
plt.show()
plt.close()
plt.figure('Spring')
nx.draw_spring(G, labels=gender)
plt.title('Spring')
plt.show()
plt.close()
plt.figure('Kamada Kawai')
nx.draw_kamada_kawai(G,labels=gender)
plt.title('Kamada Kawai')
plt.show()
plt.close()
plt.figure('Circular')
nx.draw_circular(G,labels=gender)
plt.title('Circular')
plt.show()
plt.close()
plt.figure('Spectral')
nx.draw_spectral(G,labels=gender)
plt.title('Spectral')
plt.show()
plt.close()
#Medio que en todos menos el circular y un poco menos en el spectral se ven dos grupos/clusters/comunidades. Por un lado
#uno de los clusters está compuesto principalmente por delfines de género masculino (con pocos femeninos) y
#el otro presenta la situación opuesta.
adj_real=nx.to_pandas_adjacency(G)
h=0
for i in adj_real.index:
    for j in adj_real.columns:
        if i!=j and gender[i]==gender[j] and adj_real[j][i]!=0:
            h=h+1
h_real=h/2*G.number_of_edges()
#15264.0 (real)
#asignacion aleatoria de género:
distribucion=[]
for n in np.arange(0,1000):
    copia_G=lista_G.copy()#hago una copia para no modificar el real
    random.shuffle(copia_G)
    gender_rand = dict(zip(lista_N, copia_G))
    h=0
    for i in adj_real.index:
        for j in adj_real.columns:
            if i!=j and gender_rand[i]==gender_rand[j] and adj_real[j][i]!=0:
                h=h+1
    h=h/2*G.number_of_edges()
    distribucion.append(h)

#al final no voy a fitear con la gaussiana, lo dejo como comentado por si lo queremos agregar
#sino lo borramos en la version final
#mu, std = norm.fit(distribucion)#voy a fitear con una gaussiana
plt.figure('Histograma genero random (1000)')
plt.hist(distribucion,density=True,facecolor='blue', alpha=0.5, ec='black')#normalizado
#xmin, xmax = plt.xlim()
#x = np.linspace(xmin, xmax, 100)
#p = norm.pdf(x, mu, std)
#plt.plot(x, p, 'k', linewidth=1)
plt.title('Histograma homofilia genero random (1000)')
#A priori diría que sí hay homofilia.
#Como estimación visual del valor medio, yo diría que está entre 10.500 y 12.000, podríamos decir 11.000+-1000
#Para verificarlo pordemos hacer
valor_medio=np.mean(distribucion)
std=np.std(distribucion)
#no sé muy bien como calcular el p-valor pero creo que es contar la cantidad de veces
#que en la distribucion random te dio mayor a la real y dividirlo por el total:
count=0
for i in distribucion:
    if i>h_real:
        count=count+1
p=count/len(distribucion)#p=0.0, es decir, no hay ningun caso en que la homofilia shuffleada supere la real

#por lo que vimos en clase, podemos basarnos en el grado de los nodos: (en realidad lo vimos con overlap)
#para una manera gradual, podríamos eliminar los enlaces con mayor grado (si no funciona probamos con overlap):
#probé eso y no funcionó así que después probe sacando los de menor grado primero, que para mi tiene más sentido
#y funca bárbaro pero da 46 pasos (elimina 46 nodos)

#G_copia=G.copy()
#grados=sorted(G_copia.degree, key=lambda x: x[1], reverse=False)
#pasos=0
#for i,j in grados:
#    if nx.is_connected(G_copia):
#        G_copia.remove_node(i)
#        pasos=pasos+1
#gender_2={}
#for i in G_copia.nodes:
#    gender_2[i]=gender[i]
#nx.draw(G_copia,labels=gender_2)
#comparo sacando nodos de forma random

#ahora quiero probar teniendo en cuenta el genero tambien (funciona lindo también) y elimina 36 nodos
#es decir que lo hace con menor cantidad de pasos. Yo iría con este
G_copia=G.copy()
grados=sorted(G_copia.degree, key=lambda x: x[1], reverse=False)
adj_copia=nx.to_pandas_adjacency(G_copia)
ind_col=[]
for i,j in grados:
    ind_col.append(i)
df_orden=pd.DataFrame(index=ind_col,columns=ind_col)
for i in adj_copia.index:
    for j in adj_copia.columns:
        df_orden[j][i]=adj_copia[j][i]
for i in df_orden.index:
    for j in df_orden.columns:
        if i!=j and gender[i]==gender[j] and df_orden[j][i]!=0 and nx.is_connected(G_copia):
            df_orden.drop(i,inplace=True)
            df_orden.drop(columns=j, inplace=True)
            G_copia.remove_node(i)
            break
gender_2={}
for i in G_copia.nodes:
    gender_2[i]=gender[i] 
plt.figure()
nx.draw(G_copia,labels=gender_2)
tamanio=[len(c) for c in sorted(nx.connected_components(G_copia), key=len, reverse=True)]#[17, 9]
pasos=G.number_of_nodes()-G_copia.number_of_nodes()#36


#de forma aleatoria:
tamanio_r=[]
pasos_r=[]
for n in np.arange(0,1000):
    print(n)
    G_r=G.copy()
    nodos=list(G_r.nodes()).copy()
    random.shuffle(nodos)
    count=0
    for i in nodos:
        if nx.is_connected(G_r):
            G_r.remove_node(i)
            count=count+1
    pasos_r.append(count)
    tamanio_r.append([len(c) for c in sorted(nx.connected_components(G_r), key=len, reverse=True)])

#primero podemos descartar aquellos casos que quedarons con más de 3 componentes
#Después podemos ver todos los casos que quedaron que tienen tamaños "comparables" y ver la cantidad de pasos
#que se dieron. 
mayor_2=0
steps=[]
for i in np.arange(0,len(tamanio_r)):
    if len(tamanio_r[i])>2:
        mayor_2=mayor_2+1
    elif tamanio_r[i][0]>2 and tamanio_r[i][1]>2:
        steps.append(pasos_r[i])
#grupos de tamaño 2 con 1 grupo de 1 o 2 nodos:
total_1=len(tamanio_r)-len(steps)-mayor_2 #740 (75% aprox)
#237 aprox con longitud mayor que 2
plt.figure('Pasos')
plt.hist(steps,density=False,facecolor='blue', alpha=0.5, ec='black')
plt.title('Pasos')
prom=np.mean(steps)#18 pasos aprox.


#%%

################################################################################
#                               PUNTO 3 
################################################################################

path = 'C:/Users/Mati/Documents/GitHub/Redes/tc1 data/' #Colocar la ruta donde se encuentra alojado el archivo

file = 'as-22july06.gml'


red_internet = nx.read_gml(path + file) #Abrimos el grafo de networkx desde el archivo .gml

'''
Inciso a)
'''

grados_list = [] #Inicializamos una lista que contendrá los grados de cada nodo en la red

nodos = red_internet.nodes()

for nodo in nodos:
    
    grado = red_internet.degree(nodo)
    
    grados_list.append(grado) #Agregamos el grado de cada nodo a la lista de grados

grados = np.asarray(grados_list)

print('El máximo grado alcanzado por un nodo dentro de la red es ', max(grados))  

#%%

#Bins lineales:
    
bins_lineal = np.arange(max(grados_list))

#Hacemos el histograma completo, con todos los posibles grados para los nodos de la red:
hist, bins_completo = np.histogram(grados, bins = bins_lineal, normed = True) #Con normed = True representa probabilidad

bins = bins_completo[:-1]

plt.bar(x = bins, height = hist) 

plt.xlabel('grado', fontsize = 16)

plt.ylabel('Probabilidad(grado)', fontsize = 16)

plt.show()

#%%
'''
Se observa que la mayor densidad de probabilidad se encuentra en una región de bins mucho menor al máximo.
Dado que el máximo grado que alcanzó un nodo dentro de la red es tan alto como raro de ver, decidimos truncar
el histograma. Para ello, se calculó el área debajo de la curva de probabilidades sumando el histograma paso a paso,
con tamaño de bin 1. Haciendo esto, partiendo desde 1 bin y luego agregando bins de forma iterativa, se buscó el bin
a partir del cual el área debajo de la curva es del 99%, indicando que el resto de los bins contienen simplemente
outlayers que no representan el grado típico de los nodos de la red.
'''

contador = 1

while np.sum(hist[:contador]) < 0.99: #Pedimos que pare de sumar cuando el área bajo la curva supere 0.99
    
    print(contador, np.sum(hist[:contador]))
    
    contador += 1

max_bin = contador #Este será el bin hasta el cual realizaremos el histograma.

bins_recorte = np.arange(max_bin)

hist_recorte, bins_completo_recorte = np.histogram(grados, bins = bins_recorte, normed = True)

bins_recorte = bins_completo_recorte[:-1]

plt.grid(axis = 'y')

plt.bar(x = bins_recorte, height = hist_recorte, tick_label = bins_recorte) 

plt.xlabel('grado', fontsize = 16)

plt.ylabel('Probabilidad(grado)', fontsize = 16)

plt.show()

'''
Ahora si es posible visualizar mejor la forma de la curva, luego de haber retirado los outlayers
'''
#%% 
#Bins logarítmicos:

'''
Otra manera de visualizar los datos, teniendo en cuenta la distribución irregular de los bins, es utilizar una escala
logarítmica para los mismos.
'''    

bins_log = np.log10(np.arange(1, max(grados_list)))

grados_log = np.log10(grados)

#Hacemos el histograma completo, con todos los posibles grados para los nodos de la red:
hist_log, bins_log_completo = np.histogram(grados_log, bins = bins_log+1)

#calculamos el área bajo la curva de distribución para normalizar:
area = 0

for i in range(len(hist_log)):
    
    area_barra = hist_log[i] * bins_log_completo[i]
    
    area += area_barra
    
hist_log_norm = hist_log / area #normalizamos

bins_log = bins_log_completo[:-1]

plt.grid(axis = 'y')

plt.bar(x = bins_log, height = hist_log_norm, width = np.diff(bins_log_completo))

plt.xlabel('Log(grado)', fontsize = 16)

plt.ylabel('Probabilidad(Log(grado))', fontsize = 16)

plt.show()


'''
Vemos ahora que no es necesario truncar el histograma para visualizarlo correctamente. Poniéndolo en escala logarítmica,
las diferencias abruptas se achican, dando un panorama completo de la distribución.
'''


#%%

'''
Inciso b)
'''

'''
A partir del paper: M.E.J. NEWMAN, Power laws, Pareto distributions and Zipf’s law, Contemporary Physics, Vol. 46, 
No. 5, September–October 2005, 323 – 351
se propone una ley de potencias:
'''
  
def power(x, xmin, alpha):
    
    '''
    xmin es el valor a partir del cual vale la ley de potencias
    '''
    
    C = (alpha-1)* xmin**(alpha-1)
    
    return C * x**(-alpha) 

hist_no_norm, bins_no_norm = np.histogram(grados, bins = bins_recorte) #Obtenemos el histograma recortado sin normalizar.

alpha = ig.power_law_fit(hist_no_norm, xmin = 1, method = 'discrete', return_alpha_only = True)

xmin = 1 #Sale de ver la línea anterior con return_alpha_only = False

x = np.linspace(1, max_bin, 1000)

power_fit = power(x, xmin, alpha)


#Comparamos el histograma con el fiteo:

plt.bar(x = bins_recorte, height = hist_recorte, tick_label = bins_recorte)

plt.plot(x, power_fit, color = 'red', label = 'Ley de Potencias') 

plt.xlabel('grado', fontsize = 16)

plt.ylabel('Probabilidad(grado)', fontsize = 16)

plt.legend(fontsize = 15)

plt.show()


#%%

################################################################################
#                               PUNTO 4 
################################################################################









