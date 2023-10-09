# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 10:35:35 2023

@author: Arthur Bonamigo

Dynamic Systems: Exercise 1

Double Pendulum

Arthur Bonamigo - NUSP: 12549391
João Müller ----- NUSP: 10705975
Leonardo Groot -- NUSP: 12624482
Pedro Marino ---- NUSP: 12553952
"""
import sympy as sp
from sympy import *
from sympy.physics.mechanics import *
import sympy.physics.mechanics.system as system
from sympy.physics.vector import *
init_vprinting(use_unicode = True)

#coordenadas e velocidades generalizadas:
q1, q2 = dynamicsymbols('theta1 theta2')
q1d, q2d = dynamicsymbols('theta1 theta2', 1) 
q1dd, q2dd = dynamicsymbols('theta1 theta2', 2)

#Definindo o espaço:
N = ReferenceFrame('N') #Absoluto
A1= N.orientnew('A1', 'Axis', (-(pi/2 -q1), N.z))
A2 = N.orientnew('A2', 'Axis', (-(pi/2-q2), N.z))
n1 = N.x * N.dcm(A1) #versor barra 1
n2 = N.x * N.dcm(A2) #versor barra 2
A1.set_ang_vel(N, q1d * N.z)
A2.set_ang_vel(N, q2d * N.z)

#Propriedades e parâmetros:
m1, m2, g, l1, l2 = symbols('m1 m2 g l1 l2') 

#Pontos de interesse:
O = Point('O') 
P = O.locatenew('P', l1 * A1.x)
G1= O.locatenew('G1', (l1/2) * A1.x) #Baricentro 1
Q = P.locatenew('Q', l2 * A2.x)
G2= P.locatenew('G2', l2/2 * A2.x) #Baricentro 2
O.set_vel(N, 0*N.z)

#Inércia:
I1 = (inertia(A1, m1/3*l1**2, m1/3*l1**2, m1/3*l1**2), O)
I2 = (inertia(A2, m2/3*l2**2, m2/3*l2**2, m2/3*l2**2), P)

#Barras:
B1 = RigidBody("barra_1", G1, A1, m1, I1)
B2 = RigidBody("barra_2", G2, A2, m2, I2)

#Energia Cinética:
T = simplify(kinetic_energy(N, B1, B2))

#Energia potencial:
B1.potential_energy = m1*g*(l1/2)*(1-cos(q1))
B2.potential_energy = m2*g*((l2/2)*(1-cos(q2)) + l1*(1-cos(q1)))
V = potential_energy(B1, B2)

#Método de Lagrange:
L = simplify(Lagrangian(N, B1,B2))
LM = LagrangesMethod(L, [q1, q2])
L_eq1 = simplify(LM.form_lagranges_equations())[0,0]
L_eq2 = simplify(LM.form_lagranges_equations())[1,0]
M = LM.mass_matrix

#Equações substituídas:
lamb, wp = symbols('lambda omega_p')

L_eq1_format_0 = (L_eq1.subs(l2, l1 / lamb)).subs(m2, m1 / lamb)
L_eq2_format_0 = (L_eq2.subs(l2, l1 / lamb)).subs(m2, m1 / lamb)

L_eq1_format = (L_eq1_format_0).subs(l1, g/wp**2)
L_eq2_format = (L_eq2_format_0).subs(l1, g/wp**2)

#Equações finais simplificadas:
L1 = simplify(factor(L_eq1_format))*(6*lamb**2*wp**4/g**2/m1)
L2 = simplify(factor(L_eq2_format))*(6*lamb**3*wp**4/g**2/m1)

#Sistema Linear em q1dd e q2dd:
solucao = factor(simplify(solve([L1, L2], [q1dd, q2dd])))
q1dd = solucao[q1dd]
q2dd = solucao[q2dd]


#Sistema:
F2 = Matrix(4, 1, [q1d, q2d, LM.forcing[0], LM.forcing[1]])
bodies = [B1, B2]
mass_matrix = comb_implicit_matrix = Matrix([[0, 0, 1, 0], [0, 0, 0, 1],
[ M[0,0], M[0,1] , 0, 0], [M[1,0], M[1,1], 0, 0]])
F1 = mass_matrix.inv() * F2
states = (q1, q2, q1d, q2d)
coord_idxs = (0,1)
spped_idxs = (2,3)
alg_con = [1, 2]

system1 = system.SymbolicSystem(states, F1, alg_con = alg_con,
                                bodies = bodies)

system2 = system.SymbolicSystem(states, F2, mass_matrix=comb_implicit_matrix,
                                alg_con= alg_con,
                                coord_idxs = coord_idxs,
                                bodies = bodies)


#Linearização Automática:
linearizer = LM.to_linearizer(q_ind =[q1,q2], qd_ind = [q1d, q2d])
A, B = linearizer.linearize(op_point = O, A_and_B=True, simplify=True)

#Linearização Manual:

L1_lin=((factor(L1.subs(sin(q1), q1))).subs(cos(q1 - q2), 1)).subs(3*sin(q1 - q2)*q2d**2, 0)
L2_lin = ((factor(L2.subs(sin(q2), q2))).subs(cos(q1 - q2), 1)).subs(3*lamb*sin(q1 - q2)*q1d**2, 0)

#Expressões desacopladas:
eq1 = factor(L1_lin - 3/2*L2_lin, q1dd)/lamb
eq2 = simplify(L1_lin-(2/3)*(lamb + 3)*L2_lin)

M_lin = Matrix(2,2, [2*lamb**2+6*lamb, 3, 3*lamb, 2])
K = lamb*(wp**2)*Matrix(2,2, [3*lamb + 6, 0, 0, 3])
w, w1, w2 = symbols('omega omega1 omega2')
Ç = M_lin - w**2*K
expr = collect(Ç.det(), w**2)
autovalores = solve(expr, w**2)
w1 = sqrt((autovalores[0].subs(lamb, 1.5)).subs(wp, 1))
w2 = sqrt((autovalores[1].subs(lamb, 1.5)).subs(wp, 1))





