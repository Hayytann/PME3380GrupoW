# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 01:07:16 2023

@author: Arthur Bonamigo
"""

import numpy as np
import matplotlib.pyplot as plt

# Input Parameters

v_wind = 14 # m/s
r = 63 # m
rho = 1.2 # kg/m^3

# Vectors
angular_speeds = np.linspace(0.1, 5, 100000) # rad/s
lambdas = angular_speeds*r/v_wind
betas = np.array([0 , 2,  5 , 10, 12.5, 15, 20]) # degrees
P_nom = 5e+06 # W

def power_coefficient(beta_i, lambda_k):
  factor = 1/(lambda_k + 0.08*beta_i)-0.035/(beta_i**3+1)
  Cp_i = 0.22*(116*factor - 0.4*beta_i-5)*np.exp(-12.5*factor)
  return Cp_i

Cp_matrix = np.zeros((len(betas), len(lambdas)))
for i in range(len(betas)):
  for k in range(len(lambdas)):
    Cp_matrix[i,k] = power_coefficient(betas[i], lambdas[k])


# PLOT 1: Power Coefficient Curves

Figure, axes = plt.subplots()
for i in range(len(betas)):
  label = f"beta = {betas[i]}"
  plt.plot(angular_speeds, Cp_matrix[i,:], label = label)
axes.set_ylim(0, 0.45)
plt.legend(loc = 'upper right')
plt.grid()
# axes.set_title("Coeficiente de Potência para Ângulos de Passo")
axes.set_xlabel("Velocidade de rotação do rotor (rad/s)")
axes.set_ylabel("Coeficiente de Potência (Cp)")

# Optimal Point for Wind Turbine
list_beta_12 = Cp_matrix[0].tolist()
cp_max = np.max(list_beta_12)
omega_opt = angular_speeds[list_beta_12.index(cp_max)] #rad/s
freq_opt = omega_opt/(2*np.pi)*60 # rpm
lambda_opt = omega_opt*r/v_wind

# Cp for Nominal Power Generation:
cp_nom = round(P_nom/(0.5*np.pi*r**2*v_wind**3), 3)
list_beta = np.round(Cp_matrix[0].tolist(), 3).tolist()
omegas = np.round(angular_speeds, 3)
omega_nom = angular_speeds[list_beta.index(cp_nom)]

# PLOT 2: Operation Point Displacement from Optimal
lambda_op = 4.5578
Figure2, ax = plt.subplots()
plt.plot(lambdas, Cp_matrix[0])
ax.set_title("Desvio da operação da turbina em relação ao ponto ótimo")
ax.set_xlabel("Tip speed ratio")
ax.set_ylabel("Coeficiente de Potência (Cp)")
ax.set_ylim(0, 0.45)
ax.set_xlim(0, 15)
cp_op = power_coefficient(0, lambda_op)
plt.plot(lambda_op, cp_op, "o", 
          label = f"Ponto de operação: Cp = {round(cp_op, 3)}")
plt.plot(lambda_opt, cp_max, "o",
          label = f"Ponto ótimo: Cp = {round(cp_max, 2)}")
plt.legend(loc = "lower right")
plt.grid()

# Outputs for given wind speed: 
print("velocidade do vento: ", v_wind)
print("Cp máximo: ", cp_max)
print("omega ótimo: ", omega_opt)
print("frequencia em rpm: ", freq_opt)
print("lambda ótimo: ", round(lambda_opt, 4))