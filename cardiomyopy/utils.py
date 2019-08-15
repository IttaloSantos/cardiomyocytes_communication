#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from numpy import linalg
from scipy.integrate import odeint
from cardiomyopy import probability_graphics as pg
import csv

def F(x, V1, V2, V_treshold, E1, E2, alpha1, alpha2):
    y = np.zeros(2)

    voltage = ((V2-V_treshold)/(V1-V_treshold))
    E1_coef = (E1**2 - 1)/E1
    E2_coef = (E2**2 - 1)/E2

    y[0] = (-2 - 2*voltage)*x[0]**2 + ( E1_coef*(E2**2 + 1)/((1-E2**2)*np.sqrt(alpha2)) + (2+voltage)*(E1**2 + 1)/E1 )*x[0] + ( E1_coef*(-2*E2)/((1-E2**2)*np.sqrt(alpha2)) )*x[0]*x[1] - 2

    y[1] = ( (E2**2 + 1)/(E2*np.sqrt(alpha2)) - voltage*E2_coef*((E2**2 + 1)/(1-E1**2)) )*x[1] + ( voltage/(1-E1**2)*E2_coef )*x[0]*x[1] - 2/np.sqrt(alpha2)

    return y

def JF(x, V1, V2, V_treshold, E1, E2, alpha1, alpha2):

    y = np.zeros((2,2))

    voltage = ((V2-V_treshold)/(V1-V_treshold))
    E1_coef = (E1**2 - 1)/E1
    E2_coef = (E2**2 - 1)/E2

    y[0,0] = (-2)*(2 + 2*voltage)*x[0] + ( E1_coef * ((E2**2 + 1)/((1-E2**2)*np.sqrt(alpha2))) + (2+voltage)*((E1**2 + 1)/E1) ) + ( E1_coef*(-2*E2)/((1-E2**2)*np.sqrt(alpha2)) )*x[1]
    y[0,1] = ( E1_coef*(-2*E2)/((1-E2**2)*np.sqrt(alpha2)) )*x[0]
    y[1,0] = ( voltage/(1-E1**2)*E2_coef )*x[1]
    y[1,1] = ( (E2**2 + 1)/(E2*np.sqrt(alpha2)) - voltage*E2_coef*((E2**2 + 1)/(1-E1**2)) ) + ( voltage/(1-E1**2)*E2_coef )*x[0]

    return y

def non_linear_equations_solve(V1, V2, V_treshold, E1, E2, alpha1, alpha2):
    # Calcula a solução do sistema não-linear aplicando o método de Newton

    x = np.ones(2)*0.0000001

    for i in range(100):
        delta = -np.linalg.inv(JF(x, V1, V2, V_treshold, E1, E2, alpha1, alpha2)).dot(F(x, V1, V2, V_treshold, E1, E2, alpha1, alpha2))
        x = x + delta

    return x[0], x[1]

def combination(N):
    # Cria uma matriz com todas as combinações possíveis de três números cuja soma seja igual a N

    for p in range(1, 4, 1):
        if N < 2 or p == 1:
            combinations = np.array([N, 0, 0])
            combinations = np.vstack([combinations, [0, N, 0]])
            combinations = np.vstack([combinations, [0, 0, N]])

        elif p == 2:
            for i in range(1, N, 1):
                combinations = np.vstack([combinations, [N-i, i, 0]])
                combinations = np.vstack([combinations, [N-i, 0, i]])
                combinations = np.vstack([combinations, [0, N-i, i]])

        elif p == 3:
            for j in range(1, N-1, 1):
                for i in range(1, N-j, 1):
                    combinations = np.vstack([combinations, [N-i-j, i, j]])

        else:
            raise Exception("Numero de combinacao invalido")

    return combinations

def probability_ODE(p, t, b1, a1, b2, a2):
    # Resolve equações de diferenciais do primeiro grau

    dpdt1 = b1*(1 - (p[0] + p[1])) - a1*p[0]
    dpdt2 = b2*(1 - (p[0] + p[1])) - a2*p[1]

    dpdt = [dpdt1,dpdt2]

    return dpdt

def differential_equation_solve(p0, Vj, beta1_vj, alpha1_vj, beta2_vj, alpha2_vj):
    # Resolve equações de diferenciais do primeiro grau

    return odeint(probability_ODE, p0, Vj, args=(beta1_vj, alpha1_vj, beta2_vj, alpha2_vj)) # Resolve a equação diferencial

def dump_results(Pf):
    # Salva os resultados em um arquivo CSV e gera o gráfico

    with open('/home/ittalo/Documentos/Projeto Cardiomyocytes/Codigos/data_Pf.csv', mode='w') as csv_file:
        csv_writer = csv.writer(csv_file, delimiter=',', quotechar=' ', quoting=csv.QUOTE_MINIMAL)
        for i in Pf:
            csv_writer.writerow(i)

    csv_file.close()
    pg.probability_graphic('/home/ittalo/Documentos/Projeto Cardiomyocytes/Codigos/data_Pf.csv')
