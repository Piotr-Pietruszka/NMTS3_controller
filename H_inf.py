# NMTS 3 - H_inf controller

import numpy as np
from matplotlib import pyplot as plt
import scipy.linalg as linalg
import control


def find_gamma(starting_gamma=5.0):
    gamma = starting_gamma
    stop_cond = True
    while stop_cond is True:
        A_x = A1 - D12 * np.dot(B2, C1)
        eigv_x, B_x = linalg.eigh(np.dot(B2, np.transpose(B2)) - np.power(gamma, -2) * np.dot(B1, np.transpose(B1)))
        R_x = linalg.inv(np.identity(B_x.shape[0]) * (eigv_x + 0.1))
        Q_x = np.real(linalg.sqrtm((1 - D12 * D12) * np.dot(np.transpose(C1), C1)))
        Q_x = np.dot(Q_x, np.transpose(Q_x))
        X_inf = linalg.solve_continuous_are(A_x, B_x, Q_x, R_x)

        A_y = A1 - D21 * np.dot(B1, C2)
        eigv_y, B_y = linalg.eigh(np.dot(np.transpose(C2), C2) - np.power(gamma, -2) * np.dot(np.transpose(C1), C1))
        R_y = linalg.inv(np.identity(B_y.shape[0]) * (eigv_y + 1))
        Q_y = np.real(linalg.sqrtm((1 - D21 * D21) * np.dot(B1, np.transpose(B1))))
        Q_y = np.dot(Q_y, np.transpose(Q_y))
        Y_inf = linalg.solve_continuous_are(A_y, B_y, Q_y, R_y)

        eigv, _ = linalg.eigh(np.dot(X_inf, Y_inf))
        spectral_radius = np.max(np.abs(np.real(eigv)))
        gamma -= 0.1
        if spectral_radius > np.square(gamma):
            stop_cond = False
    return gamma

# 171842,172118
a = 2
b = 8
variant_a = True

# Plant
# a)
if variant_a:
    z_0 = -abs(b - a)
# b)
else:
    z_0 = 0.5 * (b + a)
A = 0.02
# plant
# G = sig.lti([1, -z_0], [1, a+b, a*b])
s = control.TransferFunction.s
G = control.TransferFunction([1, -z_0], [1, a + b, a * b])
print(G)
print(G.pole())

# Bode plot of plant tf
mag_G, phase_G, w_G = control.bode(G, dB=True, deg=True, Plot=True)
plt.title("G bode plots")
plt.show()

# Open loop without control - Sensitivity
L = G
S = (1 + G) ** (-1)
print(S)
mag_S, phase_S, w_S = control.bode(S, dB=True, deg=True, Plot=True)
plt.title("S bode plots")
plt.show()

# GCP
print(f"z_0: {z_0}")  # w3db < 0.5x_0, gdy z_0 - zero nieminimalnofazowe
w3db = 0.7
Ms = 1.2

Wp = (s / Ms + w3db) / (s + A * w3db)

P_gcp = np.array([[Wp, -Wp * G],
                  [1, -G]])
plt.show()

# H_inf
A1 = np.asarray([[0.0, 1.0, 0.0],
                 [-a * b, -(a + b), 0.0],
                 [z_0, -1.0, A * w3db]])
B1 = np.asarray([[0.0], [0.0], [1.0]])
B2 = np.asarray([[0.0], [1.0], [0.0]])
C1 = np.asarray([[z_0, -1.0, 1.0]])
C2 = np.asarray([[z_0, -1.0, 0.0]])
D11 = 0.0
D12 = 0.0
D21 = 1.0
D22 = 0.0

gamma = find_gamma()
print(gamma)

# Closed loop - no controller
T = 1 - S
control.bode(T, dB=True, deg=True, Plot=True)
plt.title("T bode plots")
plt.show()








