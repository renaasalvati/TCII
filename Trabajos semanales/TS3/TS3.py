#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 14 19:36:48 2023

@author: renata
"""

# Inicialización e importación de módulos

# Módulos externos
import math
import numpy as np
import scipy.signal as sig

from pytc2.sistemas_lineales import analyze_sys

#Datos de enunciado
a_max = 1
a_min = 12

fp = 1500
fs = 3000

nf = 2*math.pi*fp #norma de frecuencia

wf = (2*math.pi*fp)/nf
ws = (2*math.pi*fs)/nf

e2 = 10**(a_max/10) - 1

for aux_n in range(2,5):
    aux_a_min = 10*np.log10(1 + e2*ws**(2*aux_n))
    
    if aux_a_min > a_min:
        n = aux_n
        break

#Verificación con funciones de bajo nivel
aux_r = np.roots([-1, 0, 0, 0, 0, 0, 1/e2]) #busco las raíces de T(s)*T(-s)
r = aux_r[np.real(aux_r) < 0] #me quedo con las que se encuentran en el semiplano izquierdo

num, den = sig.zpk2tf([], r, np.sqrt(1/e2))
    
# #Verificación con funciones de alto nivel
# nf2 = e2**(-1/2/n) #segunda norma de frecuencia (e⁻(1/n))

# z, p, k = sig.buttap(n)
# num, den = sig.zpk2tf(z, p, k)
# num_n, den_n = sig.lp2lp(num, den, e2**(-1/2/n)) #polinomios renormalizados con wb

analyze_sys(sig.TransferFunction(num, den), ['MP'])