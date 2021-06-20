# NMTS 3 - H_inf controller

import numpy as np
from numpy.linalg import inv
from matplotlib import pyplot as plt
import scipy.signal as sig
import control











# 171842,172118
a = 2
b = 8

# Plant
# a)
z_0 = abs(b-a)
# b)
z_0 = 0.5*(b+a)
A = 0.02
# plant
# G = sig.lti([1, -z_0], [1, a+b, a*b])
G = control.TransferFunction([1, -z_0], [1, a+b, a*b])
print(G)
print(G.pole())


# Bode plot of plant tf
mag_G, phase_G, w_G = control.bode(G, dB=True, deg=True, Plot=True)
plt.title("G bode plots")

# Bode magnitude plot
# plt.figure()
# plt.semilogx(w_G, mag_G)
# plt.xlabel('w[rad/s]')
# plt.ylabel('dB')
# plt.title("G magnitude plot")
# Bode phase plot
# plt.figure()
# plt.semilogx(w_G, phase_G)
# plt.xlabel('w[rad/s]')
# plt.ylabel('deg')
# plt.title("G phase plot")

plt.show()


# Open loop without control - Sensitivity
L = G
S = 1 + G
print(S)
mag_S, phase_S, w_S = control.bode(S, dB=True, deg=True, Plot=False)
plt.figure()
plt.semilogx(w_S, mag_S)
plt.xlabel('w[rad/s]')
plt.ylabel('dB')
plt.title("S amplitude characteristics")
plt.show()

