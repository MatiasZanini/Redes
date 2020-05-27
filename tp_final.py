# -*- coding: utf-8 -*-
"""
Created on Wed May 27 00:36:30 2020

@author: Mati
"""

import requests
from bs4 import BeautifulSoup as bs
import numpy as np
import matplotlib.pyplot as plt


#%%

path = r'C:\Users\Mati\Desktop\prueba_1.html'

a = bs(open(path), 'html')

soup = bs(a.content,'html')

f = soup.find('div', attrs={'class': '_4-u3 _5sqi _5sqk'})
likes=f.find('span',attrs={'class':'_52id _50f5 _50f7'}) #finding span tag inside class