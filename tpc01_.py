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


# Rutas a los archivos con los datos de las redes
drive_yeast_Y2H_txt = r'D:\nuestras carpetas\Mati\Libros y apuntes\Redes\Codigos\TPC 1\tc1 data\yeast_Y2H.txt'
drive_yeast_AP_MS_txt = r'D:\nuestras carpetas\Mati\Libros y apuntes\Redes\Codigos\TPC 1\tc1 data\yeast_AP-MS.txt'
drive_yeast_LIT_txt = r'D:\nuestras carpetas\Mati\Libros y apuntes\Redes\Codigos\TPC 1\tc1 data\yeast_LIT.txt'


#Inicializamos los datos:  
lista_de_enlaces_1 = abrir_txt(drive_yeast_Y2H_txt)
lista_de_enlaces_2 = abrir_txt(drive_yeast_AP_MS_txt)
lista_de_enlaces_3 = abrir_txt(drive_yeast_LIT_txt)

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

--------------------------------------------------------

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

