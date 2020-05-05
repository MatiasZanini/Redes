# -*- coding: utf-8 -*-
"""
                **Introducción a Redes Complejas en Biología de Sistemas**
                        Trabajo Computacional 1 (Entrega 05/05)


Grupo: Camila Sanz, Matías Zanini, Debora Copa
"""
################################################################################
#                                 PAQUETES 
################################################################################

import pandas as pd       #DB
import numpy as np        #matemática, simil maltab 
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import igraph as ig
import random
import math
from sklearn import linear_model
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

'''
En los diagramas de venn podemos ver que hay poca coherencia entre las redes. Podemos atribuir la falta
de coherencia a la naturaleza de las mismas. Por un lado en la red Y2H se encuentran enlazadas 
aquellas proteínas que interactuan entre sí, sin embargo, esta interacción puede no darse en el organismo 
del cual provienen. Además estos experimentos se estudian en el núcleo de la levadura 
y muchas proteínas no llegan al mismo. Por otro lado el método para generar la red AP-MS presenta 
el inconveniente de los complejos proteícos (proteínas que pueden adosarse a otras pero que no presentan 
relación funcional con las mismas). Por último la red de proteínas LIT, basada en las referencias 
bibliográficas puede presentar muchas contradicciones en sí misma, por ejemplo, un paper podría contener 
la información de que la proteína X e Y son opuestas y apreceran enlazadas en el grafo. 
Mientras que otro paper puede contener la información de X-Z cumplen la misma función y nuevamente 
apareceran enlazadas.
La red Y2H es la que presenta mayor cobertura con respecto a la cantidad de proteínas mientras que la
red AP-MS presenta mayor cobertura respecto a los enlaces.
'''

#%%

################################################################################
#                               PUNTO 2 
################################################################################

'''
Inciso a)
'''

path = 'D:/Redes 2020/TC01_data/' #Colocar la ruta donde están guardados los archivos

G = nx.read_gml( path + 'dolphins.gml')#Cargamos los nodos y enlaces

with open('D:/Redes 2020/TC01_data/dolphinsGender.txt') as f:#Cargamos los géneros
    
    gender=f.readlines()

for i in np.arange(len(gender)):
    
    gender[i]=gender[i].replace('\n','').split('\t')

gender=dict(gender)#transformamos en un diccionario

plt.figure('Original')
nx.draw(G,labels=gender)
plt.title('Original')
plt.show()
plt.close()

plt.figure('Fruchterman Reingold')
nx.draw(G, labels=gender,pos=nx.fruchterman_reingold_layout(G))
plt.title('Fruchterman Reingold')
plt.show()
plt.close()

plt.figure('Spring')
nx.draw(G, labels=gender,pos=nx.spring_layout(G))
plt.title('Spring')
plt.show()
plt.close()

plt.figure('Circular')
nx.draw(G, labels=gender,pos=nx.circular_layout(G))
plt.title('Circular')
plt.show()
plt.close()

plt.figure('Spectral')
nx.draw(G, labels=gender,pos=nx.spectral_layout(G))
plt.title('Spectral')
plt.show()
plt.close()

'''
En la mayoría de los layout testeados (exceptuando el circular y en menor grado el spectral), se distinguen dos
grupo o clusters. Uno se encuentra compuesto principalemente por delfines de género masculino, mientras que 
el otro presenta una heterogeneidad mayor entre los géneros, pero prevalecen los delfines femeninos. 
También pueden verse enlaces entre ambos grupos, es decir que las comunidades están conectadas.
En base a estas observaciones creemos que el layout "Spring" provee una visualización de la red que permite
distinguir entre los dos grupos y la conexiones entre ellos.
En este algoritmo, los enlaces actúan como resortes generando una fuerza atractiva entre nodos y 
los nodos actuán como objetos que se repelen produciendo una fuerza repulsiva, este juego
entre fuerzas continua hasta que los nodos alcazan porsiciones cercanas a las de equilibrio.
'''
#%%

'''
Inciso b)
'''
adj_real=nx.to_pandas_adjacency(G) #Matriz de adjacencia

h=0

for i in adj_real.index: #Calculamos la homofilia de la red, ignoramos el grupo NA del cual no sabemos el género
    
    for j in adj_real.columns:
        
        if i!=j and gender[i]!='NA' and gender[i]==gender[j] and adj_real[j][i]!=0:
            h=h+1

h_real=h/(2*G.number_of_edges()) #0.6037735849056604
print('Homofilia de la red: '+str(h_real))
#Asignamos aleatoriamente el género a los nodos (manteniendo la cantidad original de f,m y NA)

lista_N=[]#lista con el nombre de cada nodo

lista_G=[]#lista con el género de cada nodo

for i in gender:
    
    lista_N.append(i)
    
    lista_G.append(gender[i])

distribucion=[] #Lista de valores de homofilia para la asignación aleatoria de género

for n in np.arange(0,1000):
    
    copia_G=lista_G.copy()#Hacemos una copia de la lista de géneros para no modificar la real
    
    random.shuffle(copia_G)#Shuffleamos los géneros de forma aleatoria
    
    gender_rand = dict(zip(lista_N, copia_G)) #Generamos de nuevo el diccionario, esta vez con los generos shuffleados
    
    h=0
    
    for i in adj_real.index:
        
        for j in adj_real.columns:
            
            if i!=j and gender_rand[i]!='NA' and gender_rand[i]==gender_rand[j] and adj_real[j][i]!=0:
                h=h+1
                
    h=h/(2*G.number_of_edges())
    
    distribucion.append(h)

#Ploteamos los valores de homofilia que obtuvimos para cada asignación aleatoria

plt.figure('Histograma genero random (1000)')
plt.hist(distribucion,density=False,facecolor='blue', alpha=0.5, ec='black')
plt.title('Histograma homofilia genero random (1000)')
plt.show()
plt.close()

'''
Observando la distribución de esta variable en la asignación aleatoria de género se puede ver que la mayoría
de los valores se encuentran entre 0.42 y 0.47, es decir, podemos estimar que cuando no existe un vínculo entre
la topología de la red y el género: h=0.44+-0.03.
Debido a que h=0.6 (en la red real), podemos decir que esta red presenta homofilia.
'''
#Calculamos el valor medio esperado y el desvío standar para verificar o correjir la estimación:
valor_medio=np.mean(distribucion)
print('Valor medio: '+str(valor_medio))

std=np.std(distribucion)
print('Error: '+str(std))

#Calculamos el p-valor considerando la cantidad de veces que la homofilia dio mayor 
#en el grafo shuffleado que en el real:

count=0

for i in distribucion:
    
    if i>=h_real:
        count=count+1

p=count/len(distribucion)
print('P valor: '+str(p))

'''
Es posible que 1000 sean pocas iteraciones, probamos hacerlo con 10.000 y conseguimos p=0.0003. Con lo cual
estimamos que p=O(10^-4).
'''
#%%

'''
Inciso c)

Probamos eliminando aquellos nodos cuyo grado es menor teniendo en cuanta el género de cada nodo. Es decir, 
como posteriormente vimos que la red presentaba homofilia queremos lograr separla en dos comunidades, una
compuesta mayoritariamente por delfines de género masculino y otro mayoritariamanete por femeninos.
'''

G_copia=G.copy() #Hacemos una copia del grafo para no eliminar nodos del original

grados=sorted(G_copia.degree, key=lambda x: x[1], reverse=False)#ordenamos los grados de menor a mayor

for i,j in grados:
    
    vecinos=[n for n in G_copia.neighbors(i)]
    
    if gender[i]!='NA' and nx.is_connected(G_copia):#salteamos los que son de género desconocido
        
        for n in vecinos:
            
            if gender[i]!=gender[n] and G_copia.degree(i)<=G_copia.degree(n):#eliminamos si son de diferent género y el grado es menor
                G_copia.remove_node(i)
                break

gender_2={}

for i in G_copia.nodes:
    
    gender_2[i]=gender[i] 

plt.figure()
nx.draw(G_copia,labels=gender_2)
plt.show()
plt.close()

tamanio=[len(c) for c in sorted(nx.connected_components(G_copia), key=len, reverse=True)]#[19, 16]
print('Tamaño de las componentes: '+str(tamanio))

pasos=G.number_of_nodes()-G_copia.number_of_nodes()#27
print('Cantidad de pasos: '+str(pasos))

#Eliminando nodos de forma aleatoria:

tamanio_r=[]

pasos_r=[]

for n in np.arange(0,1000):
    
    G_r=G.copy()#de nuevo, hacemos una copia para no modificar el original
    
    nodos=list(G_r.nodes()).copy()
    
    random.shuffle(nodos)
    
    count=0
    
    for i in nodos:
        
        if nx.is_connected(G_r):
            G_r.remove_node(i)
            count=count+1
    
    pasos_r.append(count)
    
    tamanio_r.append([len(c) for c in sorted(nx.connected_components(G_r), key=len, reverse=True)])

#Primero podemos descartar aquellos casos que quedaron con más de 2 componentes

#Segundo podemos eliminar todos aquellos grupos que quedaron con órdenes de magnitud diferentes

mayor_2=0

pasos_r2=[]

tamanio_r2=[]

for i in np.arange(0,len(tamanio_r)):
    
    if len(tamanio_r[i])>2:
        mayor_2=mayor_2+1
    
    elif math.floor(math.log(tamanio_r[i][0], 10))==math.floor(math.log(tamanio_r[i][1], 10)):
        pasos_r2.append(pasos_r[i])
        tamanio_r2.append(tamanio_r[i])

#Grupos con componentes de diferentes órdenes de magnitud:

total_1=len(tamanio_r)-len(pasos_r2)-mayor_2 #(75% aprox)

print('Resultados de la eliminación aleatoria de nodos:')
print('Iteraciones con componentes de diferentes órdenes de magnitud: '+str(total_1))
print('Iteraciones con más de 2 componentes: ' + str(mayor_2))
print('Promedio de pasos realizados para las iteraciones con 2 componentes de igual orden: '+str(np.mean(pasos_r2)))# 16 aprox

diff=[]
for i in tamanio_r2:
    diff.append(abs(i[0]-i[1]))

print('Diferencia promedio entre las componentes del mismo orden:' +str(np.mean(diff)))# 16 aprox

'''
De las 1000 interaciones realizadas (eliminando nodos de manera aleatoria) quedaron aproximadamente 20 (2%) 
de los casos en los cuales la red se separó en componentes de tamaños del mismo orden. Si bien la red se separó
en una cantidad de pasos menor que la no aleatoria (16 en promedio contra 27 del original) se puede ver que
la diferencia entre los tamaños de las componentes también es superior (16 en promedio contra 3 del original).
Tomando el p-valor como el inciso anterior p=0.02 (<0.05). En base a estos resultados, nos parece que la estrategia
de eliminar secuencialmente los nodos de menor grado priorizando en eliminar aquellos nodos enlazados con
nodos de diferente género otorga resultados mejores la eliminación azarosa aunque nos hubiese gustado encontrar
una manera de separar la red con menor cantidad de pasos.
'''


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

'''
Una primera aproximación para obtener una medida de esta propiedad es obteniendo la fracción de enlaces de la red cuyos
extremos pertenecen a la misma categoría. Sin embargo esta medida es en cierta forma incompleta, ya que para casos triviales
(sólo nodos de un mismo tipo) obtendremos un valor alto. Por lo cual es conveniente sustraer de este valor la asortatividad
que se esperaría para una red independiente de la categoría.

Muchas veces nos interesan características escalares de los nodos (con valores ordinales), como por ejemplo el grado. 
El grado de un nodo  k_i , es una medida que en sí misma nos da información de la importancia de los nodos y de la 
estructura de la red. Por ejemplo pensando en la centralidad de un nodo, para el caso de las citas que recibe un artículo 
científico, es más importante cuanto más citado (conectado) está. Sin embargo, en algunos casos puede resultar de interés 
considerar el grado de centralidad de los vecinos de este nodo, ya que si estos también son centrales, éste resulta de 
mayor importancia. Esto da lugar a características interesantes de la distribución de los nodos en la red.

Modelo Barabasi
Una forma de estudiar la asortatividad es mediante un gráfico que muestre la tendencia entre el grado promedio de los nodos 
vecinos  k_(nn)(k) a los nodos de grado  k  en función de la secuencia de los grados. Si la tendencia de esta curva puede 
aproximarse con una ley de potencia se puede describir muy bien con la ecuación de Barabasi que indica:

k_nn[k] = a.k^μ 

O bien en su forma logarítmica (la cual permite hacer un ajuste lineal ya que la ley de potencia tiene una tendencia 
exponencial decreciente)

log(k_(nn) [k]) = log(a) + μ.log(k) 

Modelo Newman
Otra forma de obtener la asortatividad de grado de los nodos de la red es mediante un coeficiente, que deberá obtenerse en 
forma distinta a la asortatividad por categorías. En este tipo de casos intentar obtener el valor de la asortatividad en 
forma similar al caso de características sin orden (agrupando en bins por rangos de grado), nos llevaría al error de 
considerar características totalmente diferentes (o totalmente iguales) entre grupos, cuando en realidad no lo son y 
perderíamos así también la cercanía de los elementos entre grupos. Es por ello que se prefiere usar la covarianza 
cov (k_i, k_j)  como medida representativa de enlaces que unen a los nodos  i,j  donde  k_i  y  k_j  son variables 
aleatorias que representan el valor de el grado de cada nodo. Si normalizamos respecto a la máxima covarianza (es decir 
cuando  k_i = k_j ), obtendremos la fórmula definida por Newman del coeficiente de asortatividad  r  (equivalente a la 
correlación de Pearson).

r = ∑_(ij) (A_(ij) − k_i k_j/2m)k_i k_j / ( ∑_(ij) (k_i δ_(ij) − k_i k_j/2m)k_i k_j ) 

En este caso obtendremos valores positivos de correlación cuando los pares de nodos de la red tengan en general grados 
parecidos (asortativo) y valores negativos cuando los nodos se unen a otros con grados muy diferentes (disasortatividad) 
Esta medida de asortatividad nos da cuenta también de las características estructurales de la red. Por ejemplo, se esperan 
redes asortativas en las redes sociales, donde la gente se relaciona con sus parecidos y forman grupos. Debido a esto se 
generan distribuciones de alto grado muy concentradas en un núcleo y rodeado de una periferia poco densa.

Para casos de redes disasortativas se generan enlaces entre nodos de grados muy diferentes creando estructuras del tipo 
estrella, con una estructura más uniforme a lo largo de toda la red.
'''
#%%
# FUNCIONES

################################################################################
#    Función definida al inicio del script
################################################################################
def abrir_txt(lista_de_enlaces_txt):
    archivo=open(lista_de_enlaces_txt)
    data=[]
    for linea in archivo:
        linea=linea.strip()
        columna=linea.split()
        data.append(columna)    
    return data

################################################################################
#     obtener tabla de grados y matriz de adyacencias para archivos txt
################################################################################
def gettxt(txtfile):

  # leemos la red txt de la ruta y archivo esepcificado en datafile
  lista_de_enlaces_ = abrir_txt(txtfile)

  # Tabla de Grados 
  Red = nx.Graph()
  Red.add_edges_from(lista_de_enlaces_) #lista_de_enlaces_ es la que obtuvimos aplicando la función abrir_txt a alguno de los .txt. En networkx, no es necesario agregar primero los nodos y luego los enlaces. Podemos pasar los enlaces y agrega los nodos automáticamente.

  grados = Red.degree()                   # devuelve grado de cada nodo, como diccionario
  df_k = pd.DataFrame.from_dict(grados)   # convertimos diccionario en data frame
  df_k.columns = ['Nodos', 'Grado']       # definimos nombres de columnas
  
  # Matriz de adyacencias 
  matriz_adyacencia = nx.to_pandas_adjacency(Red) # devuelve matriz de adyacencia como dataframe de la libreria pandas
  A = matriz_adyacencia.to_numpy()  # la convertimos en matriz de numpy por si queremos hacer cuentas
  # retorna la red especificado en la ruta y archivo esepcificado en datafile
  return df_k, A

################################################################################
#      obtener tabla de grados y matriz de adyacencias para archivos gml
################################################################################
def getgml(filename):

  # extraemos el archivo con la red
  Red = nx.read_gml(filename)
      
  # Tabla de Grados 
  grados = Red.degree()                   # devuelve grado de cada nodo, como diccionario
  df_k = pd.DataFrame.from_dict(grados)   # convertimos diccionario en data frame
  df_k.columns = ['Nodos', 'Grado']       # definimos nombres de columnas

  # Matriz de adyacencias 
  matriz_adyacencia = nx.to_pandas_adjacency(Red) # devuelve matriz de adyacencia como dataframe de la libreria pandas
  A = matriz_adyacencia.to_numpy()                # la convertimos en matriz de numpy por si queremos hacer cuentas

  # retorna la red GML especificado en la ruta y archivo esepcificado en datafile
  return df_k, A


################################################################################
#                obtner grados promedio de vecinos para cada grado
################################################################################
def getknn(df_k, A):

  Knn = np.zeros(len(df_k))   # inicialización de grados promeio de vecinos

  # iteramos para cargar el grado promedio de vecinos a cada nodo i
  for i in range(0,len(df_k)):
    nn_i = A[i,:]==1                        # vector boolean de nodos conectados a i,
    Knn[i] = df_k[nn_i]["Grado"].mean()     # grado promedio de los vecinos al nodo i

  df_k['Knn promedio'] = Knn                # Nueva columna
  Knn_prom = df_k.groupby(['Grado']).mean() # promedio de los Knn para mismo grado 
  Knn_prom.index.name = 'Grado'
  Knn_prom.reset_index(inplace=True)
  
  return Knn_prom


################################################################################
#    Regrsión lineal para estimar los parametros 'a' y 'mu' del modelo Barabasai
################################################################################
def knn_barabasai(Knn_prom):

  # definimos el grado y Knn logaritmicamente para obtener una relación lineal 
  # Knn = a*K^mu   --- log -->   log(Knn) = log(a) + mu * log(K)
  logk = np.log(Knn_prom['Grado'].values)
  logknn = np.log(Knn_prom['Knn promedio'].values)
  M = np.vstack([logk, np.ones(len(logk))]).T

  #Solves the equation A x = b by computing a vector x that minimizes the squared Euclidean 2-norm || b - A x ||^2_2
  mu, loga = np.linalg.lstsq(M, logknn, rcond=None)[0]   #a = exp(loga)

  return mu, loga, logk, logknn


################################################################################
#                              Modelo Newman
################################################################################
# calculo de coeficiente de newman

def r_newman(df_k, A):

  # vector de grados para distintas potencias
  k = np.asarray(df_k["Grado"], dtype=np.float )
  k2 = k**2
  k3 = k**3

  # Sumas parciales
  S1 = np.sum(k)
  S2 = np.sum(k2)
  S3 = np.sum(k3)
  Se = (k.dot(A)).dot(k)

  # coeficiente de newman
  r = ( S1*Se - S2**2 ) / ( S1*S3 - S2**2)

  return r


#ruta donde se encuentran los archivos descargados:
path='D:/Redes 2020/TC01_data/'

papers = path+'netscience.gml'# netscience_gml

internet = path+'as-22july06.gml'# as22july06_gml

proteinas1 = path+'yeast_Y2H.txt'# yeast_Y2H_txt

proteinas2 = path+'yeast_AP-MS.txt'# yeast_AP_MS_txt

#%%

'''
Por practicidad, este ejercicio lo hicimos iterando sobre todas las redes e imprimiendo los
resultados en la consola. Los análisis y comentarios se encuentran al final del script.
'''
# iterar sobre todas las redes de ejemplo
for n in range(0,4):

  # Seleccionar red
  if n==0: 
    redsel = 'Internet'
    redfile = 'as-22july06.gml'
    item = 'a'
    filename = internet 
    df_k, A = getgml(filename) 
  elif n==1:
    redsel = 'Papers'
    redfile = 'netscience.gml'
    item = 'a'
    filename = papers 
    df_k, A = getgml(filename) 
  elif n==2:
    redsel = 'Proteinas1'
    redfile = 'yeast_Y2H.txt'
    item = 'b'
    filename = proteinas1 
    df_k, A = gettxt(filename)  

  elif n==3:
    redsel = 'Proteinas2'
    redfile = 'yeast_AP-MS.txt'
    item = 'b'
    filename = proteinas2
    df_k, A = gettxt(filename)  

  #Inciso i)
  # obtenemos la tabla de grados y matriz de adyacencias de la red
  #df_k, A = getnet(filename)  

  print('******************************************************************************')
  print(' Ejercicio 4 (' + item + '): Red de ' + redsel + ' (' + redfile + ')' )
  print('******************************************************************************')
  # obtenemos los grados de vecinos promedio
  Knn_prom = getknn(df_k, A) 
  Knn_prom = Knn_prom.drop([0]) # eliminamos los grados 0

  #Inciso ii)
  # graficamos los grados de vecinos promedio
  
  print(redsel + ' a-ii) Grado promedio de vecindad para nodos de cierto grado (Knn)')
  Knn_prom.plot(x='Grado', y='Knn promedio', style='o')  
  plt.title('Grado vs Knn promedio')  
  plt.xlabel('Grado')  
  plt.ylabel('Knn promedio')  
  plt.show()

  #Inciso iii)
  #obtenemos los parametros de Knn mediante el estimador de Barabasai
  
  mu, loga, logk, logknn = knn_barabasai(Knn_prom)
  a = np.exp(loga) # parametro 'a' en escala lineal

  # curva del estimador de barabasai
  Knn_barab = loga + mu * logk 

  # grafico de los puntos obtenidos de la red y el ajuste lineal
  print(redsel +' a-iii) Parámetros Barabasai estimados:\na = ' + str(a) + '\nmu = ' + str(mu))

  plt.plot(logk, logknn, 'o', label='Original data', markersize=5)
  plt.plot(logk, Knn_barab, 'r', label='Fitted line')
  plt.legend()
  plt.title('Grado vs Knn promedio')  
  plt.xlabel('log(Grado)')  
  plt.ylabel('log(Knn_promedio)')  
  plt.show()
  
  #Inciso iv)
  #parametros estimados

  r = r_newman(df_k, A)
  print(redsel + ' a-iv) Parámetro de newman (r): ' + str(r) + '\n')

#%%
'''
Gráfico de la red de colaboraciones científicas 
'''

Red = nx.read_gml(papers)
plt.figure() 
nx.draw_kamada_kawai(Red, 

        width = .5, # Ancho de los enlaces
        node_color = 'blue', # Color de los nodos
        edge_color = 'violet', # Color de los enlaces
        node_size = 5, # Tamaño de los nodos
        with_labels = False # True o False, Pone o no el nombre de los nodos
       )
plt.title('Red de colaboraciones científicas')

#%%
'''
En términos generales, se puede verificar que existe una relación entre la tendencia de la curva  Knn(k)  
y el coeficiente de Newman  r , notando que éste resulta positivo para tendencias positivas de la curva, 
mostrando la asortatividad, mientras que para tendencias negativas el coeficiente de Newman resulta también negativo, 
indicando disarsotatividad en ese caso.

Red de Internet:
Observamos que la curva Knn(k) posee un decaimiento exponencial, lo cual indica una tendencia negativa 
y por lo tanto la disasortatividad.
La mayoría de los nodos tienen un grado bajo y unos pocos nodos en la cola de la curva de grado alto que se conectan 
sólo con nodos de grado bajo. También se puede observar la misma curva pero en escala logarítmica (ambos ejes) y 
la recta cuya ordenada al origen es  log(a) (a=575.7) y pendiente mu=−0.44 ,
Por otro lado el cálculo del coeficiente  r  arroja un valor negativo de  r=−0.198 , 
corroborandose la correlación negativa entre los grados de los nodos de la red. 
En este caso se esperaría una estructura del grafo de la red de internet en forma de estrella.

Red de colaboraciones científicas:
Para el caso de la red de colaboraciones científicas, en la curva  Knn(k) se puede apreciar que existen nodos de grado 
bajo que se relacionan en general con sus parecidos, lo cual indica la asortatividad, sin embargo luego 
del grado k=10 esta relación deja de valer y los nodos de la red se conectan con nodos de grado diferente 
aunque con una leve tendencia creciente pero mucho menos asortativa. 
En este ejemplo se obtuvieron como parámetros de ajuste lineal (que no es el mejor modelo de ajuste en este caso) 
con  a=3.55  y  mu=0.30 . Luego, el coeficiente de Newman resulta también positivo r=0.46. 
Por lo tanto se pueden encontrar redes asortativas, no solo en los casos de redes sociales, 
mostrando que las publicaciones tienden a citar artículos de su mismo campo (están agrupados).

Red de proteinas (yeast_Y2H.txt)
La curva  Knn(k)  muestra una tendencia negativa, por lo tanto disasortativa, aunque muy dispersa. 
Según este gráfico, podemos pensar en que esta dispersión muestra que muchos nodos de grado alto se conectan con 
nodos de grado alto y bajo y viceversa. Por lo cual la ley de potencias no da suficiente información en este caso 
más que una idea de una leve disasortatividad, cuyos parámetros del ajuste lineal resultan a=17,58 y mu=-0,2. 
Podemos sin embargo observar que la baja correlación de la asortatividad de la red se manifiesta con un coeficiente 
de Newman muy bajo, arrojando un valor de r=-0.05.
Respecto a la estructura de la red puede verse que los nodos están más conectados no formándose el conglomerado central 
tan marcado como en el otro caso.

Red de proteinas (yeast_AP_MS.txt)
La curva  Knn(k)  muestra una tendencia marcadamente positiva, mostrando que es asortativa, ya que se observan 
nodos de grado bajo que se conectan con nodos de grado alto por lo menos hasta el grado 70 aproximadamente, 
luego esa tendencia baja levemente, pero se mantiene positiva si se considera todo el conjunto. 
En el ajuste lineal sobre esta curva, que también se muestra en escala logarítmica, se obtienen a=4.47 y mu=0.6. 
Por el lado del coeficiente de Newman, se puede apreciar una alta correlación positiva r=0.6, 
verificando su comportamiento asortativo. En este caso al observar la estructura se puede ver claramente 
un grupo de alta densidad cerca del centro

'''
#%%