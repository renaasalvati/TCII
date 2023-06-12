#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  4 16:09:10 2023

@author: renata
"""

# Inicialización e importación de módulos

# Módulos externos
import numpy as np
import scipy.signal as sig
import matplotlib as mpl
import math

from matplotlib import pyplot as plt
from pytc2.sistemas_lineales import analyze_sys, pretty_print_lti, bodePlot, pzmap, tf2sos_analog, pretty_print_SOS


fig_sz_x = 13
fig_sz_y = 7
fig_dpi = 80 # dpi
fig_font_size = 13

mpl.rcParams['figure.figsize'] = (fig_sz_x, fig_sz_y)
mpl.rcParams['figure.dpi'] = fig_dpi
plt.rcParams.update({'font.size':fig_font_size})

#Datos de enunciado para filtro pasa altos
a_max = 3
a_min = 20

fp1 = 1600*(10**3)
fp2 = 2500*(10**3)
fs1 = 1250*(10**3)
fs2 = 3200*(10**3) 

w0 = 2*np.pi*np.sqrt(fp1*fp2)

nf = w0 #norma de frecuencia

#normalización de frecuencias angulares
w0_n = w0/nf
wp1_n = (2*np.pi*fp1)/nf
wp2_n = (2*np.pi*fp2)/nf
ws1_n = (2*np.pi*fs1)/nf
ws2_n = (2*np.pi*fs2)/nf

Bw = wp2_n - wp1_n #ancho de banda
Q = w0_n/Bw #factor de selectividad

#conversión de parámetros para un filtro pasa bajos equivalente
Wp1_n = round(Q*(wp1_n**2 - w0_n**2)/wp1_n, 3)
Wp2_n = round(Q*(wp2_n**2 - w0_n**2)/wp2_n, 3)
Ws1_n = round(Q*(ws1_n**2 - w0_n**2)/ws1_n, 3)
Ws2_n = round(Q*(ws2_n**2 - w0_n**2)/ws2_n, 3)

Wp_n = math.fabs(Wp1_n)

#verifico para Ws1 y Ws2 cuál es la más chica en módulo y elijo esa para trabajar con el prototipo
if math.fabs(Ws1_n) < math.fabs(Ws2_n):
    Ws_n = math.fabs(Ws1_n)
else:
    Ws_n = math.fabs(Ws2_n)     

print(f"Frecuencias angulares normalizadas del pasa-bajos prototipo: \nwp = {Wp_n}, ws = {Ws_n}\n")    
    
#Obtención de epsilon y n
e2 = round(10**(a_max/10) - 1, 1)

for aux_n in range(2,5):
    aux_a_min = 10*np.log10(1 + e2*Ws_n**(2*aux_n))
    
    if aux_a_min > a_min:
        n = aux_n
        break
        
print(f"Parámetros del filtro pasa-bajos prototipo: \ne² = {e2}, n = {n}")

#transferencias de cada etapa del pasa-bajos prototipo 
num1_lp_n = num2_lp_n = [1.]
den1_lp_n = [1., 1.]
den2_lp_n = [1., 1., 1.]

#numeradores y denominadores resultantes de la transformación de cada etapa del pasa-bajos
num1_bp_n, den1_bp_n = sig.lp2bp(num1_lp_n, den1_lp_n, w0_n, Bw)
num2_bp_n, den2_bp_n = sig.lp2bp(num2_lp_n, den2_lp_n, w0_n, Bw)

sos1_bp = tf2sos_analog(num1_bp_n, den1_bp_n)
sos2_bp = tf2sos_analog(num1_bp_n, den2_bp_n)

num_bp_n = np.polymul(num1_bp_n, num2_bp_n)*3.1623
den_bp_n = np.polymul(den1_bp_n, den2_bp_n)

sos_bp = tf2sos_analog(num_bp_n, den_bp_n)

print("Primera etapa del pasa-bajos al pasa-banda")
pretty_print_lti(num1_bp_n, den1_bp_n)
pretty_print_SOS(sos1_bp, mode='omegayq') #factorizando

print("\nSegunda etapa del pasa-bajos al pasa-banda")
pretty_print_lti(num2_bp_n, den2_bp_n)
pretty_print_SOS(sos2_bp, mode='omegayq') #factorizando

print("\nArmando la transferencia completa")
pretty_print_lti(num_bp_n, den_bp_n)
pretty_print_SOS(sos_bp, mode='omegayq') #factorizando

#-----
num_lp_n = [1.]
den_lp_n = [1., 2., 2., 1.]

#transferencias completas del pasa-bajos y pasa-banda
H_lp_n = sig.TransferFunction(num_lp_n, den_lp_n)
H_bp_n = sig.TransferFunction(num_bp_n, den_bp_n)

sos_bp = tf2sos_analog(num_bp_n, den_bp_n)

analyze_sys(sos_bp)

