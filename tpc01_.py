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

# Evitar acumulación de mensajes de warning en el display
import warnings  
warnings.filterwarnings("ignore")


#%%

################################################################################
#                               PUNTO 1 
################################################################################

def abrir_txt(lista_de_enlaces_txt):

    '''
    Abre los datos de las redes guardados en .txt
    '''
    
    archivo = open(lista_de_enlaces_txt)
    
    data=[]
    
    for linea in archivo:
    
        linea=linea.strip()
        
        columna=linea.split()
        
        data.append(columna)    
    
    return data

# Red = nx.read_gml(drive_as_22july06_gml)

path='D:/nuestras carpetas/Mati/Libros y apuntes/Redes/Codigos/TPC 1/tc1 data/'
filename=['yeast_AP-MS','yeast_LIT','yeast_Y2H']
df=pd.DataFrame(index=filename,columns=['N','L','<k>','k_max','k_min','rho','<C_i>','C','D'])
ne=pd.DataFrame(index=filename,columns=['nodos','enlaces'])#lo vamos a usar después para los diagramas de venn
for file in filename:
    with open(path+file+'.txt') as f:#abrimos el archivo
        data=f.readlines()
    
    for i in np.arange(len(data)):#transformamos en una lista de tuples
        data[i]=tuple(data[i].replace('\n','').split('\t'))

    G=nx.Graph()#inicializamos el grafo
    G.add_edges_from(data)#agregamos los nodos y enlaces
    ne.loc[file]=[G.nodes(),G.edges()]
    plt.figure()
    plt.title(file)
    nx.draw(G, 
        width = 2, # Ancho de los enlaces
        node_color = 'lightblue', # Color de los nodos
        edge_color = 'black', # Color de los enlaces
        node_size = 10, # Tamaño de los nodos
        with_labels = False # True o False, Pone o no el nombre de los nodos
       )
    #plt.savefig('D:/Redes 2020/TC01_ejercicios/1a_'+file+'.png')
    plt.show()
    plt.close()
    df['N'][file]=G.number_of_nodes()
    df['L'][file]=G.number_of_edges()
    df['<k>'][file]=2*G.number_of_edges()/G.number_of_nodes()
    grados=[m for n,m in G.degree()]
    df['k_max'][file]=np.max(grados)
    df['k_min'][file]=np.min(grados)
    df['rho'][file]=2*G.number_of_edges()/(G.number_of_nodes()*(G.number_of_nodes()-1))
    #creo que esto se puede hacer con nx.density()
    df['<C_i>'][file]=nx.average_clustering(G)
    df['C'][file]=nx.transitivity(G)
    if nx.is_connected(G):
        df['D'][file]=nx.diameter(G)
    else:
        df['D'][file]='infinito'

df.to_csv('D:/nuestras carpetas/Mati/Libros y apuntes/Redes/Codigos/TPC 1/tc1 data/1b.csv')

plt.figure()
nodos_venn = venn3([set(ne['nodos'][filename[0]]), set(ne['nodos'][filename[1]]),set(ne['nodos'][filename[2]])],set_labels = ('yeast_AP-MS', 'yeast_LIT', 'yeast_Y2H'))
plt.title('Diagrama de Venn (nodos)')
#plt.savefig('D:/Redes 2020/TC01_ejercicios/1d_nodos.png')

plt.figure()
nodos_venn = venn3([set(ne['enlaces'][filename[0]]), set(ne['enlaces'][filename[1]]),set(ne['enlaces'][filename[2]])],set_labels = ('yeast_AP-MS', 'yeast_LIT', 'yeast_Y2H'))
plt.title('Diagrama de Venn (enlaces)')
#plt.savefig('D:/Redes 2020/TC01_ejercicios/1d_ejes.png')

'''
algunas cosas para decir
El grafo es maso menos parecido para todas las redes, hay un conjunto de proteínas de levadura en el
medio y un conjunto afuera que lo rodea. En terminos de nodos, la última red tiene aprox 500 más que los otros.
Podemos ver en la cantidad de enlaces y el grado que la primer red, donde los enlaces se posicionan si
2 proteínas interactúan es superior a los otros (lo mismo ocurre con el grado máximo). La densidad es similar
para las 3 redes, un poco superior para la primera (coíncide con los enlaces y el grado). El coeficiente
de clustering mustra que hay mayor cantidad de nodos que cumplen la clausura transitiva en la primer red, es
decir si la proteína X interactúa con la proteína Y y la Z, es probable que las proteínas Y y Z interactúen entre
sí (de hecho como es parecido a 1/2 te diría que la mitad de los nodos cumplen la clausura transitiva. 
En cuanto a la segunda red, vemos que el coeficiente de clustering es menor (aproximadamente 1/3 o 1/4 de 
los nodos cumplen CT), suponemos uqe esto se debe a que una proteína puede cumplir más de una función
y aquellas proteínas cuyos vecinos están enlazados por diferentes funciones, no se enlazan entre sí.
Por último el coeficiente de clustering es menos en la 3 red, es decir, que no necesariamente si
las proteínas X-Y aparecen en un paper y X-Z aparece en otro, no hay una gran probabilidad de que 
Y-Z aparezcan en otro paper, principalmente porque los paper pueden estar hablando de diferentes temáticas.
Algo similar se puede deducir del coeficiente de clustering global.
Dado que las redes no son conexas, por convención el diametro es infinito.
'''



# ruta_archivos = 'D:\nuestras carpetas\Mati\Libros y apuntes\Redes\Codigos\TPC 1\tc1 data\'

# # Rutas a los archivos con los datos de las redes
# drive_yeast_Y2H_txt = (ruta_archivos + 'yeast_Y2H.txt')
# drive_yeast_AP_MS_txt = ruta_archivos + 'yeast_AP-MS.txt'
# drive_yeast_LIT_txt = ruta_archivos + 'yeast_LIT.txt'


# #Inicializamos los datos:  
# lista_de_enlaces_1 = abrir_txt(drive_yeast_Y2H_txt)
# lista_de_enlaces_2 = abrir_txt(drive_yeast_AP_MS_txt)
# lista_de_enlaces_3 = abrir_txt(drive_yeast_LIT_txt)

# from google.colab import drive
# drive.mount('/content/drive')
#%%
'''
Inciso a)
    Presentamos una comparación gráfica entre las 3 redes
'''

Red_proteinas_1 = nx.Graph()

Red_proteinas_1.add_edges_from(lista_de_enlaces_1) #lista_de_enlaces_1 es la que obtuvimos aplicando la función abrir_txt a alguno de los .txt. 

plt.figure() #abrimos una nueva figura

plt.title('Red de proteinas: yeast_Y2H')

nx.draw(Red_proteinas_1, 
#nx.draw_kamada_kawai(Red_proteinas_1, 

        width = .5, # Ancho de los enlaces

        node_color = 'blue', # Color de los nodos

        edge_color = 'violet', # Color de los enlaces

        node_size = 5, # Tamaño de los nodos

        with_labels = False # True o False, Pone o no el nombre de los nodos
       )


Red_proteinas_2 = nx.Graph()

Red_proteinas_2.add_edges_from(lista_de_enlaces_2) #lista_de_enlaces_2 es la que obtuvimos aplicando la función abrir_txt a alguno de los .txt.

plt.figure() #abrimos una nueva figura

plt.title('Red de proteinas: yeast-AP_MS')

nx.draw(Red_proteinas_2, 
#nx.draw_kamada_kawai(Red_proteinas_2, 
        width = .5, # Ancho de los enlaces

        node_color = 'violet', # Color de los nodos

        edge_color = 'blue', # Color de los enlaces

        node_size = 5, # Tamaño de los nodos

        with_labels = False # True o False, Pone o no el nombre de los nodos
       )

Red_proteinas_3 = nx.Graph()

Red_proteinas_3.add_edges_from(lista_de_enlaces_3) #lista_de_enlaces_3 es la que obtuvimos aplicando la función abrir_txt a alguno de los .txt. 

plt.figure() #abrimos una nueva figura

plt.title('Red de proteinas: yeast_LIT')

nx.draw(Red_proteinas_3, 
#nx.draw_kamada_kawai(Red_proteinas_3, 

        width = .5, # Ancho de los enlaces

        node_color = 'green', # Color de los nodos

        edge_color = 'red', # Color de los enlaces

        node_size = 5, # Tamaño de los nodos

        with_labels = False # True o False, Pone o no el nombre de los nodos
       )

plt.show()


#%%

'''
Inciso b)
    Listamos varias propiedades de las redes
'''

# i.

#Número de nodos en las 3 redes:
N_nodos_1 = Red_proteinas_1.number_of_nodes()

N_nodos_2 = Red_proteinas_2.number_of_nodes()

N_nodos_3 = Red_proteinas_3.number_of_nodes()

#ii.

#Número de enlaces en las 3 redes:
N_enlaces_1 = Red_proteinas_1.number_of_edges()

N_enlaces_2 = Red_proteinas_2.number_of_edges()

N_enlaces_3 = Red_proteinas_3.number_of_edges()


#iii.

'''
Si no se lo aclara, uno puede interpretar que la red es tanto dirigida, como no dirigida. Sin embargo, es posible
establecer un criterio basado en la matriz de adyacencia. Si la matriz de adyacencia de la red es simétrica, 
resulta intuitivo pensar que la red es no dirigida, ya que no se tiene en cuenta el sentido de los enlaces.
Vemos si la matriz es simétrica simplemente restándole su transpuesta y chequeando que se anule.
'''

matriz_adyacencia_1 = nx.to_pandas_adjacency(Red_proteinas_1) # devuelve matriz de adyacencia como dataframe de la libreria pandas

print('La matriz de Adyacencia es: ')

print(matriz_adyacencia_1) #mostramos como se ve una matriz tipica de adyacencia

matriz_adyacencia_2 = nx.to_pandas_adjacency(Red_proteinas_2)

matriz_adyacencia_3 = nx.to_pandas_adjacency(Red_proteinas_3)

Adj1 = matriz_adyacencia_1.to_numpy() # la convertimos en matriz de numpy

Adj2 = matriz_adyacencia_2.to_numpy()

Adj3 = matriz_adyacencia_3.to_numpy()

matrices = [Adj1, Adj2, Adj3]

for i in matrices:
    
    A = i - np.transpose(i)

    if np.all(A==0):
        
        print('La matriz es simétrica. Por lo tanto, consideramos la red como no dirigida')
        
    else:
        
        print('La matriz no es simétrica. La red es dirigida')
    
    





    

#%%








"""#Ejercicio 1
1) Considere las tres redes de interacción de proteínas relevadas para levadura disponibles en la
página de la materia. Se trata de: una red de interacciones binarias (yeast_Y2H.txt), de copertenencia a complejos proteicos (yeast_AP-MS.txt) y obtenida de literatura (yeast_LIT.txt)
obtenidas del Yeast Interactome Database.

a. Presente una comparación gráfica de las 3 redes.

b. Resuma en una tabla las siguientes características de dichas redes

i. El número total de nodos, N

ii. El número total de enlaces L, de la red

iii. Si se trata de una red dirigida o no-dirigida

iv. El grado medio <k> (<kin>,<kout> en caso de red dirigida), el grado máximo y
mínimo de la red.

v. La densidad de la red

vi. Los coeficientes de clustering $<Ci>$ y $C_Δ$ de la red.

vii. Diámetro de la red.


c. Teniendo en cuenta la naturaleza de las interacciones reportadas, diga si es razonable lo
que encuentra para ciertos observables calculados.

d. Construya un diagrama de Venn que permita reconocer la cobertura, especificidad y coherencia de las interacciones reportadas por los tres datasets

-----------------------------------------------------
La idea de este ejercicio es indagar en algunas de las características topológicas principales de tres redes de interacción de proteínas de la levadura de cerveza.
##Inciso (a)
En este inciso, queremos simplemente visualizar las tres redes. Para esto, primero, necesitamos generarnos las redes a partir de las listas de enlaces obtenidas de la lectura de los .txt. Luego, podemos generar las visualizaciones con networkx.
```
import networkx as nx
import matplotlib.pylab as plt # se recomienda fuertemente importar todos los paquetes en una celda aparte al principio del cuaderno

Red_proteinas_1 = nx.Graph()
Red_proteinas_1.add_edges_from(lista_de_enlaces_1) #lista_de_enlaces_1 es la que obtuvimos aplicando la función abrir_txt a alguno de los .txt. En networkx, no es necesario agregar primero los nodos y luego los enlaces. Podemos pasar los enlaces y agrega los nodos automáticamente.
nx.draw(Red_proteinas_1)
plt.show()
```
##Inciso (b)
En este inciso buscamos comparar características topológicas de las redes. El grado, coeficiente de clustering, etc. Cada una de estas características puede obtenerse con funciones propias de la librería networkx. Recomendamos que ustedes mismxs las busquen y si tienen dudas nos consultan. Una forma interesante de dar cuenta de las distintas características de las redes es mediante una tabla. Para esto, pueden usar el paquete pandas
```
import pandas as pd
# Lo más cómodo es utilizar diccionarios para cada columna
diccionario_columna_nodos = {'Redes' : [Red_proteinas_1, Red_proteinas_2, Red_proteinas_3], 'Número de nodos' : [N1,N2,N3]} # y así una llave para cada una de las características que querramos
tabla_comparativa = pd.DataFrame(data = diccionario_columna_nodos) # Existen múltiples atributos para esta función, recomendamos fuertemente incursionar en la documentación de la librería
# Para visualizar la tabla, basta con escribir el nombre de la misma en la última línea de la celda utilizada y ejecutar
tabla_comparativa
```

##Inciso (c)
En este inciso se busca discutir lo obtenido en el anterior teniendo en cuenta la naturaleza de cada una de las redes estudiadas. 

*Comentarios Debora:* 

*Las 3 redes se arman con aproximadamente los mismos nodos, pero generando enlaces en forma diferente:
La primer red pone enlaces si las proteinas interaccionan entre si.
La segunda red pone enlaces si estan ligadas a la misma funcion
La tercer red usa tecnicas de Data Mining a partir de papers relacionados con proteinas de levaduras y genera un enlace si los nombres de de diferentes proteinas están en una misma oracion. Esto trae problemas ya que por ejemplo la oracion: Proteina1 no interacciona con Proteina2, generaría un enlace en la red.*
##Inciso (d)
En este inciso se busca utilizar diagramas de Venn para estudiar qué tan distintas son las redes en el sentido de enlaces y nodos. Es decir, se sugiere realizar dos diagramas de Venn: uno para los nodos y otro para los enlaces.
"""

'''
b. Resuma en una tabla las siguientes características de dichas redes
i. El número total de nodos, N
.'''
cantidad_nodos = Red_proteinas_1.number_of_nodes()    # devuelve la cantidad de nodos
lista_nodos = list(Red_proteinas_1.nodes())           # devuelve los nombres de los nodos

print('La red tiene ' + str(cantidad_nodos) + ' nodos')
print('Los nodos son ' + str(lista_nodos))

#ii. El número total de enlaces L, de la red
cantidad_enlaces = Red_proteinas_1.number_of_edges()  # devuelve la cantidad de enlaces
lista_enlaces = list(Red_proteinas_1.edges())         # devuelve los nombres de los enlaces
print('La red tiene ' + str(cantidad_enlaces) + ' enlaces')
print('Los enlaces son ' + str(lista_enlaces))

matriz_adyacencia = nx.to_pandas_adjacency(Red_proteinas_1) # devuelve matriz de adyacencia como dataframe de la libreria pandas
print('La matriz de Adyacencia es: ')
print(matriz_adyacencia)

#iii. Si se trata de una red dirigida o no-dirigida
Adj1 = matriz_adyacencia.to_numpy() # la convertimos en matriz de numpy
A = Adj1+np.transpose(Adj1)
np.all(A==0)

#iv. El grado medio (, en caso de red dirigida), el grado máximo y mínimo de la red
#tomar la componente conexa más grande (CC max)
R_P_1_ND=Red_proteinas_1.to_undirected() #R_P_1_ND la vuelve no dirigida
largest_cc = max(nx.connected_components(R_P_1_ND), key=len) #set de nodos de la CC max
smallest_cc = min(nx.connected_components(R_P_1_ND), key=len) #set de nodos de la CC min
'''
plt.figure()
MCC=R_P_1_ND.subgraph(largest_cc) #crea una subgrafica con el conjuto largest_cc
nx.draw_kamada_kawai(MCC, 
        width = .5, # Ancho de los enlaces
        node_color = 'green', # Color de los nodos
        edge_color = 'red', # Color de los enlaces
        node_size = 5, # Tamaño de los nodos
        with_labels = False # True o False, Pone o no el nombre de los nodos
       )
plt.show()

plt.figure()
SCC=R_P_1_ND.subgraph(smallest_cc) #crea una subgrafica con el conjuto largest_cc
nx.draw_kamada_kawai(SCC, 
        width = .5, # Ancho de los enlaces
        node_color = 'green', # Color de los nodos
        edge_color = 'red', # Color de los enlaces
        node_size = 5, # Tamaño de los nodos
        with_labels = False # True o False, Pone o no el nombre de los nodos
       )
plt.show()
'''
#v. La densidad de la red
Densidad = nx.density(Red_proteinas_1)
print('La densidad de la red es: ' + str(Densidad))

#vi. Los coeficientes de clustering y CΔ de la red.
print('El coeficientes de clustering promedio es: ' + str(nx.average_clustering(Red_proteinas_1))) #https://networkx.github.io/documentation/stable/reference/algorithms/generated/networkx.algorithms.cluster.average_clustering.html#networkx.algorithms.cluster.average_clustering

#vii. Diámetro de la red
#Diametro = nx.diameter(Red_proteinas_1) #The diameter function relies on a net being strongly connected. For a weakly connected graph one could use: Diametro = nx.diameter(Red_proteinas_1.to_undirected())
# Solo me interesa la componente gigante
Diametro = nx.diameter(largest_cc)

print('El diámetro de la red es: ' + str(Diametro))


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



"""--------------------------------------------------------

La idea general de este ejercicio es, a partir de la red de sistemas autónomos de internet, adquirir la noción de distribución de grado y las distintas estrategias que pueden utilizarse para estudiar dicha distribución.
## Inciso (a)
La idea de este inciso es graficar la distribución de grado de distintas maneras para tener una apoximación a la naturaleza de dicha distribución.
Si bien existen muchas librerías que permiten generar histogramas, las dos más usuales son: numpy y matplotlib. La segunda nos permite mediante un simple comando generar el histograma y graficar al mismo tiempo. En cambio, con la primera, generamos la distribución, y después necesitamso de otra librería para graficar.


```
# Lista con los grados de los nodos (una de las posibilidades de generarla)

lista_de_grados = []
for nodo in Red.nodes()
  grado = Red.degree(nodo)
  lista_de_grados.append(grado)

# Caso 1
plt.hist(lista_de_grados) # La función plt.hist() tiene muchos atributos: cantidades de bins, normalización, etc
plt.show()

# Caso 2
hist, bins = np.histogram(lista_de_grados) # Esta función también tiene varios atributos: cantidad de bins, espaciado, etc.
# bins es un array que da cuenta de dónde comienza y dónde termina cada bin, por lo tanto, tiene un elemento más que hist, que son la cantidad de nodos que caen dentro de un determinado bin. Tenemos que definir si tomaremos el centro de los bines, el comienzo o el final. En el siguiente ejemplo tomamos el comienzo
plt.bar(x = bins[:-1], height = hist, width = np.diff(bins))
plt.show()

```
La idea es que prueben cómo se ve la distribución de grado al variar la escala con la que se grafica (lineal o logarítmica) y al variar el bineado que se utiliza, es decir: bineado lineal o bineado logaritmico (usar un tamaño de bin constante o uno que varíe).
## Inciso (b)
Se busca encontrar la relación funcional para la distribución de grado. Para esto, se recomienda utilizar la librería scipy y sus funciones para ajustar.

#Ejercicio 4
Asortatividad

a. Considere la red de colaboraciones científicas (netscience.gml) y la red de internet (asjuly06.gml). Analice si nodos de alto grado tienden a conectarse con nodos de alto grado
o por el contrario suelen conectarse a nodos de bajo grado? (i.e la red es asortativa o disortativa respecto al grado?). Para ello:

i. Determine, para nodos de grado k, cuánto vale en media el grado de sus vecinos.
[hint R: se puede estimar primero el grado medio de los vecinos de cada nodo de
la red y luego utilizar aggregate sobre esos datos, que permite estimar cantidades
sobre subconjuntos determinados de datos de acuerdo a diferentes criterios]

ii. Analizar la tendencia observada en un gráfico que consigne dicho valor knn(k)
como función del grado.

iii. Asumiendo que k_{nn}(k)=akmu
, estime el exponente de correlación a partir de
realizar una regresión de . Asegurese de graficar el fiteo en el
grafico anterior. [hint R: lm permite hacer regresiones lineales]

iv. Considere la red de colaboraciones y la de internet nuevamente Encuentre
cuantitativamente la asortatividad de la red utilizando ahora el estimador
propuesto por Newman:

Para ello tenga encuenta lo desarrollado en las eqs [8.26 – 8.29] del libro de
Newman.Como se corresponde este coeficiente con el estimado en el punto
anterior? A qué se debe?

b. Corra el script de cálculo (puntos i-iii) para las redes Y2H y AP-MS. Puede explicar lo
que observa en cuanto a la asortatividad reportada?

------------------------------------------

La idea principal de este ejercicio es estudiar si en las redes propuestas (ojo, la red de colaboraciones es una red pesada) existe asortatividad en el grado. Más en general, comprender cuál es el grado medio de los vecinos de un nodo, en función del grado de este nodo. A su vez, se pide que se repita el análisis para las redes de proteínas vistas anteriormente.

# Inciso (a)
Para este inciso, es importante entender los pasos necesarios para llevar acabo la tarea pedida. Entendemos que, en primer lugar, es recomendable trabajar con el diccionario de nodos y sus respesctivos grados antes que una lista de grados (pensar, de forma alternativa, si puede ser útil trabajar con un diccionario cuyas llaves -keys- sean los distintos grados y los valores -values- listas de nodos con determinado grado). A su vez, para tener acceso a los vecinos de un nodo, podemos hacer uso de la matriz de adyacencia de la red, pero también tenemos una función de la librería networkx que nos permite acceder a un iterable con los vecinos de determinado nodo:

 -------------------------------------------------------

En la practica, para construir knn(k) se suelen seguir los siguientes pasos. Primero, se calcula el grado
de cada nodo de la red. Luego, se los agrupa a los nodos por grado. Tercero, para cada conjunto
de nodos de determinado grado, se calcula el grado medio de los vecinos de cada nodo y luego se
lo promedia por el total de nodos que determinado grado. De esta forma, se construye la relaci ́on
knn(k).



```
# Opción 1
vecinos_nodo_i = Red.neighbors(i) # donde Red es un nx.Graph() e i un nodo cualquiera

# Opción 2
vecinos_nodo_i = Red[i] # ojo porque acá obtenemos un diccionario donde podemso tener información sobre el enlace entre el nodo i y sus vecinos (por ejemplo el peso)
```
La idea final es que estudiemos la relación entre el grado medio de los vecinos de los nodos de grado k en función de k. El estudio de esta relación se debe hacer en base a los dos modelos propuestos (Newman y Barabasai) sobre el origen de la asortatividad.

## Inciso (b)
Repetir lo anterior pero para las redes de proteínas
"""

#%%

################################################################################
#                               PUNTO 4 
################################################################################









