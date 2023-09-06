#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 16 00:19:15 2023

@author: renata
"""
# Módulos externos
import numpy as np
import scipy.signal as sig
import matplotlib as mpl

from matplotlib import pyplot as plt
from pytc2.sistemas_lineales import analyze_sys, pretty_print_lti, bodePlot, pzmap, tf2sos_analog, pretty_print_SOS

fig_sz_x = 13
fig_sz_y = 7
fig_dpi = 80 # dpi
fig_font_size = 13

mpl.rcParams['figure.figsize'] = (fig_sz_x, fig_sz_y)
mpl.rcParams['figure.dpi'] = fig_dpi
plt.rcParams.update({'font.size':fig_font_size})

#Datos de enunciado
d = 80*(10**-6)

f1 = 3*(10**3)
f2 = 2*(10**3)

nf = 1/d

w1_n = (2*np.pi*f1)/nf
w2_n = (2*np.pi*f2)/nf

print(f"Las frecuencias angulares normalizadas son:\nw1 = {w1_n} w2 = {w2_n}")

n = 3

z,p,k = sig.besselap(n, 'delay')
num, den = sig.zpk2tf(z, p, k)

H = sig.TransferFunction(num, den)
pretty_print_lti(num, den)

SOS_H = tf2sos_analog(num, den)
pretty_print_SOS(SOS_H, 'omegayq')

analyze_sys(H);