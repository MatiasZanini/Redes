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

import numpy as np        #matemática, simil maltab 
import networkx as nx
import matplotlib.pyplot as plt
import random
import math
import plfit
from collections import Counter
from scipy.optimize import curve_fit
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

path2='D:/Redes 2020/TC02/graficos/'

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
    
    plt.savefig(path2+'hist_'+filename[i]+'.png')
    
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
    
    plt.savefig(path2+'scatter_'+filename[i]+'.png')

    plt.show()
    
    plt.close()
    
#%%
'''
Inciso b)
Cualitativamente, se puede ver que la distribución de grados de la red de Internet (as-22july06) y 
de la red de proteínas (yeast_Y2H) son las que mejor se ajustan a una power law. 
Para verificar estas observaciones, podemos hacer los ajustes correspondientes (inciso d)
'''

'''
Inciso c)
Todas las redes exiben efectos de borde, por un lado, ninguna red presenta nodos con k=0. Por otro lado, 
si cualitativamente trazamos una recta lineal en los histogramas partiendo del mínimo grado en donde
podemos comenzar a trazar la recta estimamos que deberíamos ver hubs de grado superior en la red.
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
    
    plt.savefig(path2+filename[i]+'_ajuste.png')
    
    plt.show()
    
    plt.close()
    
    xmin = fit._xmin
    
    alpha = fit._alpha

    print('Red '+filename[i]+': K_min = '+str(xmin)+'; Gamma = '+str(alpha))
    
    gammas.append(alpha)
    
    kminimo.append(xmin)

#%%
'''
Inciso e)

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

#Definimos una función útil para utilizar durante el punto:
    
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
Inciso a)
'''
# Armamos la red random cuyos nodos se conectan con probabilidad constante p:

def red_erdos(p, n):
    
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
    
    return red_er






    
p = 0.2 # Probabilidad con la que se conectan dos nodos.

n = 1e4 # Cantidad de Nodos en la red

red_er = red_erdos(p, n)

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
Inciso b)
'''

def red_random(k0, n):
    
    # Creamos una nueva red del tipo cliqué de grado k0, es decir, hay k0 nodos todos enlazados entre sí:
    
    red_rand = nx.complete_graph(k0)

    for ki in range(k0, n):
        
        # Ahora, creamos una lista de nodos al azar (no repetidos) de la red con la que se conectará el nuevo nodo:
        enlaces = random.sample(list(red_rand.nodes()), k = k0)
        
        for i in enlaces:
            
            red_rand.add_edge(i, ki) # Agregamos cada enlace para el nuevo nodo.
            
    return red_rand

'''
Con esto hemos creado un algoritmo iterativo para generar una red, cuyos nodos iniciales tienen grado k0, agregando un nuevo
nodo de grado k0 uniéndose a los demás de forma aleatoria en cada paso.
'''
    
k0 = 5 # Establecer el grado inicial de la red

n = int(1e4) # Establecer el número de nodos que se desea que tenga la red.

red_rand = red_random(k0, n)


#%%

'''
Inciso c)
'''

'''
La red tipo Barabasi es similar a la generada en el item b). Sin embargo, aquí cada nodo agregado se enlazará a los demás
nodos de la red con una probabilidad que depende del grado de los mismos. Cuando mayor sea el grado, mayor la probabilidad 
de que el nuevo nodo se conecte a él. Según Barabasi: p(k_i) = k_i / sum(k_i).
'''

def red_barabasi(k0, k, n):
    
    if k>k0:
        
        raise ValueError('No se puede añadir un nodo con un grado k mayor que el inicial k0 sin repetir enlaces.')
    
    red_barab = nx.complete_graph(k0)     
    
    for ki in range(k0, n):
        
        grado_arr = np.asarray(red_barab.degree())[:,1] # Genera un array con los grados de la red
        
        probs = grado_arr / sum(grado_arr) # Array con las probabilidades de cada grado segun Barabasi: p(k_i) = ki / sum(k_i)
        
        '''
        Ahora, creamos una lista de nodos elegidos con la probabilidad dada por probs (no repetidos) de la red 
        con la que se conectará el nuevo nodo:
        '''
        enlaces = np.random.choice(np.asarray(red_barab.nodes()), size = k, p = probs, replace = False)
    
        for i in enlaces:
            
            red_barab.add_edge(i, ki) # Agregamos cada enlace para el nuevo nodo.
            
    return red_barab


# Creamos una nueva red del tipo cliqué de grado k0, es decir, hay k0 nodos todos enlazados entre sí:
    
k0 = 5 # Establecer el grado inicial de la red

n = int(1e4) # Establecer el número de nodos que se desea que tenga la red.

k = 5 # Establecer el grado de los nodos que se agregarán en cada paso. IMPORTANTE: k <= k0

red_barab = red_barabasi(k0, k, n)

# i.

'''
Ahora tenemos que comparar la cantidad m de enlaces en la red con n*k, en particular podríamos usar k=k0:
'''    

m = red_barab.number_of_edges() # Cantidad de enlaces

m_teo_b = n*k # Valor esperado para el número de enlaces en una red del tipo Barabasi

print('El número de enlaces <k> difiere del valor esperado en un', np.around(abs(m-m_teo_b)/m_teo_b * 100, 4), '%')

'''
Como vemos, ambos valores son comprables. Esto se debe a la aparición de hubs, tal como se esperaba. Los nodos con gran
cantidad de enlaces (grado alto), tienden a captar los nuevos nodos agregados a la red. Esto implica que para un número 
grande de nodos, la cantidad de enlaces se encuentre dominada por estos hubs. 

Esto quiere decir que, como en cada paso se agregó un 
nodo de grado k, y que la mayoría de los enlaces fueron a parar a dichos hubs, el tamaño de los enlaces totales en la red
resulta similar a multiplicar k por la cantidad de pasos. Además, como n>>k0>k, se tiene que, si el número de pasos es
n-k0 ---> el número de pasos será similar a n, con lo cual m será similar a k*n.
'''

#%%

'''
Inciso d)
'''


'''
Caso 1. Comparando con "as-22july06_edgelist" :
'''
# Ponemos los valores de n y m que obtuvimos para el ejercicio 1:
n1 = 22941

m1 = 48372


p = m1 / (n1*(n1-1)/2) # Calculamos el p que tendría asociado una red erdos-renyi de n nodos y m enlaces.

red_er_1 = red_erdos(p, n1) #Red erdos-renyi para estos valores de n y m

red_rand_1 = red_random(5, n1) # Elegimos k0 = 5  de forma arbitraria para inicializar la red random.

k1 = int(m1/n1) # Calculamos el grado k de cada nodo agregado a la red tipo barabasi. 

red_barab_1 = red_barabasi(5, k1, n1) # Elegimos k0 = 5 de forma arbitraria para inicializar la red random.

#%%

# Graficamos:
    
grados_er = np.asarray(red_er_1.degree())[:,1]

#grados_er_log = np.log10(grados_er)

#bins_er = np.logspace(np.min(grados_er), np.max(grados_er), 13)

bins_er = np.arange(np.min(grados_er), np.max(grados_er)+1, 1)

plt.figure()
    
hist, bins, bars = plt.hist(grados_er, bins = bins_er, align = 'mid', density = True, facecolor='blue', alpha=0.5, ec='black')

#plt.xscale('log')

plt.yscale('log')

plt.title('Distribución de grado para red tipo Erdos-Renyi simulando la red "as-22july06_edgelist" ')

plt.xlabel('k (linear scale)')

plt.ylabel('P(k) (log scale)')

plt.grid()

plt.show()


'''
Observamos una distribución del tipo Poisson cuyo máximo se alcanza en el grado k=4. Esto es esperable, ya que el grado
medio de la red es 2*m1/n1 el cual es aproximadamente 4.23.

Debido a que no aparecieron Hubs, no fue necesario emplear la escala logaritmica en el eje x (bins), ya que los grados de 
cada nodo son bajos (el máximo es 15).
'''

#%%

grados_rand = np.asarray(red_rand_1.degree())[:,1]

grados_rand_log = np.log10(grados_rand)

bins_rand = np.logspace(np.min(grados_rand_log), np.max(grados_rand_log), 20)

plt.figure()
    
hist, bins, bars = plt.hist(grados_rand, bins = bins_rand, align = 'mid', density = True, facecolor='blue', alpha=0.5, ec='black')

plt.xscale('log')

plt.yscale('log')

plt.title('Distribución de grado para red tipo Random simulando la red "as-22july06_edgelist" ')

plt.xlabel('k (log scale)')

plt.ylabel('P(k) (log scale)')

plt.grid()

plt.show()


#%%

grados_barab = np.asarray(red_barab_1.degree())[:,1]

grados_barab_log = np.log10(grados_barab)

bins_barab = np.logspace(np.min(grados_barab_log), np.max(grados_barab_log), 20)

plt.figure()
    
hist, bins, bars = plt.hist(grados_barab, bins = bins_barab, align = 'mid', density = True, facecolor='blue', alpha=0.5, ec='black')

plt.xscale('log')

plt.yscale('log')

plt.title('Distribución de grado para red tipo Barabasi simulando la red "as-22july06_edgelist" ')

plt.xlabel('k (log scale)')

plt.ylabel('P(k) (log scale)')

plt.grid()

plt.show()

'''
En este histograma, se puede observar claramente un decaimiento del tipo libre de escala. Esto era lo esperado ya que
el modelo de Barabasi propone una mejor aproximación a las redes reales que los modelos random.
'''



#######################################################################################################################
#%%

'''
Caso 2. Comparando con "netscience_edgelist" :
'''

n2 = 1450

m2 = 2727

p = m2 / (n2 *(n2 - 1)/2) # Calculamos el p que tendría asociado una red erdos-renyi de n nodos y m enlaces.

red_er_2 = red_erdos(p, n2) #Red erdos-renyi para estos valores de n y m

red_rand_2 = red_random(5, n2) # Elegimos k0 = 5  de forma arbitraria para inicializar la red random.

k2 = int(m2/n2) # Calculamos el grado k de cada nodo agregado a la red tipo barabasi. 

red_barab_2 = red_barabasi(5, k2, n2) # Elegimos k0 = 5 de forma arbitraria para inicializar la red random.

#%%

# Graficamos:
    
grados_er = np.asarray(red_er_2.degree())[:,1]

#grados_er_log = np.log10(grados_er)

#bins_er = np.logspace(np.min(grados_er), np.max(grados_er), 13)

bins_er = np.arange(np.min(grados_er), np.max(grados_er)+1, 1)

plt.figure()
    
hist, bins, bars = plt.hist(grados_er, bins = bins_er, align = 'mid', density = True, facecolor='blue', alpha=0.5, ec='black')

#plt.xscale('log')

plt.yscale('log')

plt.title('Distribución de grado para red tipo Erdos-Renyi simulando la red "netscience_edgelist" ')

plt.xlabel('k (linear scale)')

plt.ylabel('P(k) (log scale)')

plt.grid()

plt.show()


'''
Análogo al caso anterior. Se observa una curva tipo Poisson.
'''

#%%

grados_rand = np.asarray(red_rand_2.degree())[:,1]

grados_rand_log = np.log10(grados_rand)

bins_rand = np.logspace(np.min(grados_rand_log), np.max(grados_rand_log), 20)

plt.figure()
    
hist, bins, bars = plt.hist(grados_rand, bins = bins_rand, align = 'mid', density = True, facecolor='blue', alpha=0.5, ec='black')

plt.xscale('log')

plt.yscale('log')

plt.title('Distribución de grado para red tipo Random simulando la red "netscience_edgelist" ')

plt.xlabel('k (log scale)')

plt.ylabel('P(k) (log scale)')

plt.grid()

plt.show()


#%%

grados_barab = np.asarray(red_barab_2.degree())[:,1]

grados_barab_log = np.log10(grados_barab)

bins_barab = np.logspace(np.min(grados_barab_log), np.max(grados_barab_log), 12)

plt.figure()
    
hist, bins, bars = plt.hist(grados_barab, bins = bins_barab, align = 'mid', density = True, facecolor='blue', alpha=0.5, ec='black')

plt.xscale('log')

plt.yscale('log')

plt.title('Distribución de grado para red tipo Barabasi simulando la red "netscience_edgelist" ')

plt.xlabel('k (log scale)')

plt.ylabel('P(k) (log scale)')

plt.grid()

plt.show()

'''
Si bien en este histograma se encuentra, nuevamente como era esperado, un decaimiento libre de escala, se ve que los 
últimos puntos se alejan de este comportamiento. Es posible que se deba a que esta red posee muchos menos nodos que la
anterior, aumentando el error estadístico.
'''

#######################################################################################################################
#%%
'''
Caso 3. Comparando con "power_enlaces" :
'''
# Ponemos los valores de n y m que obtuvimos para el ejercicio 1:
n3 = 4941

m3 = 6594


p = m3 / (n3*(n3-1)/2) # Calculamos el p que tendría asociado una red erdos-renyi de n nodos y m enlaces.

red_er_3 = red_erdos(p, n3) #Red erdos-renyi para estos valores de n y m

red_rand_3 = red_random(5, n3) # Elegimos k0 = 5  de forma arbitraria para inicializar la red random.

k3 = int(m3/n3) # Calculamos el grado k de cada nodo agregado a la red tipo barabasi. 

red_barab_3 = red_barabasi(5, k3, n3) # Elegimos k0 = 5 de forma arbitraria para inicializar la red random.


#%%

# Graficamos:
    
grados_er = np.asarray(red_er_3.degree())[:,1]

#grados_er_log = np.log10(grados_er)

#bins_er = np.logspace(np.min(grados_er), np.max(grados_er), 13)

bins_er = np.arange(np.min(grados_er), np.max(grados_er)+1, 1)

plt.figure()
    
hist, bins, bars = plt.hist(grados_er, bins = bins_er, align = 'mid', density = True, facecolor='blue', alpha=0.5, ec='black')

#plt.xscale('log')

plt.yscale('log')

plt.title('Distribución de grado para red tipo Erdos-Renyi simulando la red "power_enlaces" ')

plt.xlabel('k (linear scale)')

plt.ylabel('P(k) (log scale)')

plt.grid()

plt.show()


'''
Análogo a los casos anteriores. Se observa una curva tipo Poisson.
'''

#%%

grados_rand = np.asarray(red_rand_3.degree())[:,1]

grados_rand_log = np.log10(grados_rand)

bins_rand = np.logspace(np.min(grados_rand_log), np.max(grados_rand_log), 20)

plt.figure()
    
hist, bins, bars = plt.hist(grados_rand, bins = bins_rand, align = 'mid', density = True, facecolor='blue', alpha=0.5, ec='black')

plt.xscale('log')

plt.yscale('log')

plt.title('Distribución de grado para red tipo Random simulando la red "power_enlaces" ')

plt.xlabel('k (log scale)')

plt.ylabel('P(k) (log scale)')

plt.grid()

plt.show()


#%%

grados_barab = np.asarray(red_barab_3.degree())[:,1]

grados_barab_log = np.log10(grados_barab)

bins_barab = np.logspace(np.min(grados_barab_log), np.max(grados_barab_log), 13)

plt.figure()
    
hist, bins, bars = plt.hist(grados_barab, bins = bins_barab, align = 'mid', density = True, facecolor='blue', alpha=0.5, ec='black')

plt.xscale('log')

plt.yscale('log')

plt.title('Distribución de grado para red tipo Barabasi simulando la red "power_enlaces" ')

plt.xlabel('k (log scale)')

plt.ylabel('P(k) (log scale)')

plt.grid()

plt.show()

'''
Nuevamente, para el caso de la red tipo Barabasi, se encuentra una distribución libre de escala, tal como se esperaba.
'''

#######################################################################################################################

#%%
'''
Caso 4. Comparando con "yeast_Y2H" :
'''
# Ponemos los valores de n y m que obtuvimos para el ejercicio 1:
n4 = 2018

m4 = 2930


p = m4 / (n4*(n4-1)/2) # Calculamos el p que tendría asociado una red erdos-renyi de n nodos y m enlaces.

red_er_4 = red_erdos(p, n4) #Red erdos-renyi para estos valores de n y m

red_rand_4 = red_random(5, n4) # Elegimos k0 = 5  de forma arbitraria para inicializar la red random.

k4 = int(m4/n4) # Calculamos el grado k de cada nodo agregado a la red tipo barabasi. 

red_barab_4 = red_barabasi(5, k4, n4) # Elegimos k0 = 5 de forma arbitraria para inicializar la red random.

#%%

# Graficamos:
    
grados_er = np.asarray(red_er_4.degree())[:,1]

#grados_er_log = np.log10(grados_er)

#bins_er = np.logspace(np.min(grados_er), np.max(grados_er), 13)

bins_er = np.arange(np.min(grados_er), np.max(grados_er)+1, 1)

plt.figure()
    
hist, bins, bars = plt.hist(grados_er, bins = bins_er, align = 'mid', density = True, facecolor='blue', alpha=0.5, ec='black')

#plt.xscale('log')

plt.yscale('log')

plt.title('Distribución de grado para red tipo Erdos-Renyi simulando la red "yeast_Y2H" ')

plt.xlabel('k (linear scale)')

plt.ylabel('P(k) (log scale)')

plt.grid()

plt.show()


'''
Análogo a los casos anteriores. Se observa una curva tipo Poisson.
'''

#%%

grados_rand = np.asarray(red_rand_4.degree())[:,1]

grados_rand_log = np.log10(grados_rand)

bins_rand = np.logspace(np.min(grados_rand_log), np.max(grados_rand_log), 18)

plt.figure()
    
hist, bins, bars = plt.hist(grados_rand, bins = bins_rand, align = 'mid', density = True, facecolor='blue', alpha=0.5, ec='black')

plt.xscale('log')

plt.yscale('log')

plt.title('Distribución de grado para red tipo Random simulando la red "yeast_Y2H" ')

plt.xlabel('k (log scale)')

plt.ylabel('P(k) (log scale)')

plt.grid()

plt.show()


#%%

grados_barab = np.asarray(red_barab_4.degree())[:,1]

grados_barab_log = np.log10(grados_barab)

bins_barab = np.logspace(np.min(grados_barab_log), np.max(grados_barab_log), 13)

plt.figure()
    
hist, bins, bars = plt.hist(grados_barab, bins = bins_barab, align = 'mid', density = True, facecolor='blue', alpha=0.5, ec='black')

plt.xscale('log')

plt.yscale('log')

plt.title('Distribución de grado para red tipo Barabasi simulando la red "yeast_Y2H" ')

plt.xlabel('k (log scale)')

plt.ylabel('P(k) (log scale)')

plt.grid()

plt.show()

'''
Si bien no sigue la tendencia libre de escala tan bien como en el caso anterior, se observa que se la red tipo Barabasi
muestra este comportamiento de manera satisfactoria.
'''

'''
Como conclusión final, se encuentra que las redes del tipo Random y Erdos-Renyi no son satisfactorias a la hora de predecir
el comportamiento de las redes reales. En ninguna de las redes estudiadas de este tipo se generaron Hubs. Sin embargo,
en las redes reales estos Hubs, o nodos de alto grado, si que aparecen. 
Por su parte, las redes de tipo Barabasi, por construcción, establecen una mayor prioridad a la conexión de nuevos nodos
con los demás de la red cuyos grados sean altos. Esto propicia no solo la aparición de Hubs, sino también un comportamiento
libre de escala, el cual se observa en las redes reales.
'''




#%%

################################################################################
#                               PUNTO 3 
################################################################################

#Usamos el mismo código que en el Punto 2, inciso c)

# Creamos una nueva red del tipo cliqué de grado k0, es decir, hay k0 nodos todos enlazados entre sí:
    
k0 = 7 # Establecer el grado inicial de la red

n = int(10000) # Establecer el número de nodos que se desea que tenga la red.

red_barab = nx.complete_graph(k0)

k = 7 # Establecer el grado de los nodos que se agregarán en cada paso. IMPORTANTE: k <= k0

t=0#paso

#Guardamos en listas los valores de los pasos y los grados para t=5 y t=95

step_5=[]

grado_5=[]

step_95=[]

grado_95=[]

if k>k0:
    
    raise ValueError('No se puede añadir un nodo con un grado k mayor que el inicial k0 sin repetir enlaces.')
    

for ki in range(k0, n):
        
    t=t+1
    
    grado_arr = np.asarray(red_barab.degree())[:,1] # Genera un array con los grados de la red
    
    probs = grado_arr / sum(grado_arr) # Array con las probabilidades de cada grado segun Barabasi: p(k_i) = ki / sum(k_i)
    
    '''
    Ahora, creamos una lista de nodos elegidos con la probabilidad dada por probs (no repetidos) de la red 
    con la que se conectará el nuevo nodo:
    '''
    enlaces = np.random.choice(np.asarray(red_barab.nodes()), size = k, p = probs, replace = False)
    
    for i in enlaces:
        
        red_barab.add_edge(i, ki) # Agregamos cada enlace para el nuevo nodo.
    
    if t>=5:
        
        grado_5.append(red_barab.degree[k0+5-1])#restamos 1 porque comenzamos a contar del nodo 0
        
        step_5.append(t)
    
    if t>=95:
        
        grado_95.append(red_barab.degree[k0+95-1])#restamos 1 porque comenzamos a contar del nodo 0
        
        step_95.append(t)

fig = plt.figure()

ax1 = fig.add_subplot(111)

ax1.scatter(step_5, grado_5, s=10, c='C0', label='t=5')

ax1.scatter(step_95,grado_95, s=10, c='C2', label='t=95')

ax1.set_xscale('log')

ax1.set_yscale('log')

plt.legend(loc='upper left')

plt.ylabel('Grado (k)')

plt.xlabel('t (paso)')

plt.savefig(path2+'kvst_5_95.png')

plt.show()

plt.close()

'''
Podemos observar que, en un principio, las curvas de la evolución temporal del grado para los nodos
seleccionados son inestables y a medida que incrementamos la cantidad de nodos en la red, adoptan un comportamiento
aproximadamente lineal (en escala logarítmica) con una pendiente similar. Esto puede observarse con
más claridad en la curva de t=5, como este nodo es agregado en pasos posteriores, se puede
observar la estabilidad de la curva para t>100.
Podemos entonces estimar el exponente (pendiente en la escala adoptada) de estas curvas.
'''
#fiteamos

#dividimos por el t0 de cada una para que las curvas comiencen del mismo x0-y0

x5=np.divide(step_5,step_5[0])

x95=np.divide(step_95,step_95[0])

def exp_func(x, a):
    
    return k * (x**a)

popt5, pcov5 = curve_fit(exp_func, x5,grado_5)

popt95, pcov95 = curve_fit(exp_func, x95, grado_95)

fig = plt.figure()

ax1 = fig.add_subplot(111)

ax1.scatter(x5, grado_5, s=10,c='C0', label='t=5')

ax1.scatter(x95,grado_95, s=10,c='C2', label='t=95')

ax1.plot(x5, exp_func(x5, *popt5), 'r-',label='Fit: a=%5.3f (t=5)' % tuple(popt5))

ax1.plot(x95, exp_func(x95, *popt95), 'r--',label='Fit: a=%5.3f (t=95)' % tuple(popt95))

ax1.set_xscale('log')

ax1.set_yscale('log')

plt.legend(loc='upper left')

plt.ylabel('Grado (k)')

plt.xlabel('t/t0 (paso)')

plt.savefig(path2+'kvst_5_95_fit.png')

plt.show()

plt.close()

'''
Como podemos ver, la pendiente de ambas curvas es similar a 0.48+-0.03 (este valor puede estar ligeramente
modificado por la forma aleatoria en que se genera la red, recomendamos incrementar la cantidad de 
de pasos para observar comportamientos más estables de las curvas, en particular de la curva
para t=95). Podemos estimar que este comportamiento, para t>t0+100 resulta independiente del nodo que tomemos.
Los hubs de la red corresponden, en general, a los nodos más viejos (aquellos agregados en los primeros
pasos, es decir, con menor t0), incrementan su conectividad en una proporción mayor que los nuevos. 
Debido a que la probabilidad de establecer un enlace entre un vértice nuevo y viejo es proporcional al grado de 
este último, se genera el efecto de "rich-get-richer", en donde aquellos nodos con mayor cantidad de conexiones 
tendrán probabilidad más alta de establecer un nuevo enlace. 
De esta forma los hubs incrementan su grado y por ende su probabilidad de conexión en cada iteración.
'''
#%%
