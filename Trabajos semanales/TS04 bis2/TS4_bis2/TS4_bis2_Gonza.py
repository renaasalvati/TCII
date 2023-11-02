import matplotlib as plt
import pytc2.sistemas_lineales as tc2
import scipy.signal as sig
import math as m

def w2Omega(w, Q):

    Omega = abs(Q*(w**2 -1)/w)

    return Omega

def ordenCheby(alpha_min, ee ,ws):

    # Otra forma
    # n = m.log(m.sqrt(4*(10**(alpha_min/10)-1)/ee))/m.log(ws + m.sqrt(ws**2 - 1)
    n = 0
    alpha_aux = 0
    while alpha_aux < alpha_min:
        n += 1
        alpha_aux = 10*m.log10(1+ee*m.cosh(n*m.acosh(ws))**2)

    return n

# Plantilla
f0          = 22e3      # Hz
fs1         = 17e3      # Hz
fs2         = 36e3      # Hz
alpha_max   = 0.5       # dB
alpha_min1  = 16        # dB
alpha_min2  = 24        # dB
Q           = 5

# Normalizo plantilla
w0  = 1
ws1 = fs1/f0
ws2 = fs2/f0

BW = w0/Q

# Armo la plantilla del pasa bajo prototipo

Omega_0  = 1
Omega_s1 = w2Omega(ws1, Q)
Omega_s2 = w2Omega(ws2, Q)

# Propongo trabajar con una transferencia Butter

# Obtengo epsilon
ee = 10**(alpha_max/10) - 1
e = m.sqrt(ee)

# obtengo el orden del filtro
n1 = ordenCheby(alpha_min1, ee, Omega_s1,)
n2 = ordenCheby(alpha_min2, ee, Omega_s2,)

n = max(n1, n2)

print(f"El orden del filtro es {n1}\n")

# Obtengo el filtro prototipo

z, p, k =sig.cheb1ap(n, alpha_max)
num, den = sig.zpk2tf(z, p, k)
tf_lp = sig.TransferFunction(num, den)

tc2.analyze_sys(tf_lp, sys_name="Filtro prototipo")

print("Transferencia del filtro pasabajos prototipo resulta:/n")
tc2.pretty_print_lti(num, den)

# factorizo la expresion
sos = tc2.tf2sos_analog(num, den)
print("\nTransferencia del prototipo factorizada en secciones de segundo orden:\n")
tc2.pretty_print_SOS(sos)


# Realizo la transformación

num, den = sig.lp2bp(num, den, w0, BW)
tf_bp = sig.TransferFunction(num, den)

tc2.analyze_sys(tf_bp, sys_name="Filtro PasaBanda")

print("Transferencia del filtro Pasabanda normalizado:")
tc2.pretty_print_lti(num, den)

# factorizo la expresion
sos = tc2.tf2sos_analog(num, den)
print("\nTransferencia del pasa banda normalizada y factorizada en secciones de segundo orden:\n")
tc2.pretty_print_SOS(sos)


