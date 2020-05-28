# -*- coding: utf-8 -*-
"""
Created on Wed May 27 00:36:30 2020

@author: Mati
"""
import os
from datetime import datetime
#import requests
from bs4 import BeautifulSoup as bs
import numpy as np
import matplotlib.pyplot as plt

#%%

'''
Funciones
'''

def nombre_post(file, h_post, n_post):
    
    h = file.find(h_post)
    
    n = file.find(n_post, h)
    
    nlen = len(n_post)
    
    n_i = n + nlen + 1 # Inicio del nombre
    
    n_f = file.find('"', n_i) # Fin del nombre

    return file[n_i : n_f]


def fecha_post(file, date_file, h_date, n_date):
    
    h = file.find(h_date)
    
    n = file.find(n_date, h)
    
    n_i = file.find('>', n) + 1 # Inicio de la fecha
    
    n_f = file.find('<', n_i) # Fin de la fecha
    
    fecha_cruda = file[n_i : n_f]
    
    if fecha_cruda.find('h') > 0:
        
        fecha_2_list = [date_file.year, date_file.month, date_file.day, date_file.hour - [int(s) for s in fecha_cruda.split() if s.isdigit()][0]]
        
        if fecha_2_list[3] < 0:
            
            fecha = [date_file.year, date_file.month, date_file.day -1 , 24 + fecha_2_list[3]]
            
        else:
            
            fecha = fecha_2_list
    
    else:
        
        fecha = fecha_cruda
    
    return fecha


def nombre_reacts(file, h_reacts, n_reacts):
    
    reacters = []
    
    nlen = len(n_reacts)
    
    h = file.find(h_reacts)
    
    h_stop = file.find('Nuevo mensaje', h)
    
    stop = h_stop - nlen -1
    
    n = file.find(n_reacts, h)
    
    while n < stop:
    
        n_i = n + nlen + 1 # Inicio del nombre
    
        n_f = file.find('"', n_i) # Fin del nombre

        nombre = file[n_i : n_f]
        
        if nombre not in ['Agregar', 'Seguir', 'Mensaje']:

            reacters.append(nombre)
            
            print(nombre)
            
        n = file.find(n_reacts, n_f + 1)

    return reacters






#%%

path = 'C:/Users/Mati/Documents/GitHub/Redes/datos/'

filename = 'fuente_completa_likes_2.html'

date_file = datetime.fromtimestamp(os.path.getctime(path+filename)) # Objeto con la fecha y hora de la obtención de los datos. 
#Atributos: .year, .month, .day, .minute, .second

'''
Tags:
'''

h_post = '<div class="sjgh65i0 l9j0dhe7 k4urcfbm du4w35lb">' # Head para encontrar al posteador

n_post = '<a aria-label=' # Head para el nombre del posteador

h_date = '<div class="qzhwtbm6 knvmm38d">' # Head para encontrar la fecha del post

n_date = '<span id="jsc_c' # Head para el valor de la fecha

h_reacts = '<div class="a8s20v7p k5wvi7nf buofh1pr pfnyh3mw l9j0dhe7 du4w35lb"><div data' # Head para encontrar los que reaccionan

n_reacts = 'aria-label=' # Head para cada nombre de los que reaccionan. Ignorar "Agregar", "Seguir", "Mensaje"
# Cortar en "Nuevo mensaje" 


file = str(bs(open(path+filename, encoding="utf8"), 'html.parser'))
    
poster = nombre_post(file, h_post, n_post)

reacters = nombre_reacts(file, h_reacts, n_reacts)

fecha = fecha_post(file, date_file, h_date, n_date)





































#%%

path = r'C:\Users\Mati\Documents\GitHub\Redes\datos\fuente_completa_likes_1.html'

a = bs(open(path, encoding="utf8"), 'html')

soup = bs(a.content,'html')

f = soup.find('div', attrs={'class': '_4-u3 _5sqi _5sqk'})
likes=f.find('span',attrs={'class':'_52id _50f5 _50f7'}) #finding span tag inside class