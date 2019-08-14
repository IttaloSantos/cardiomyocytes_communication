#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from numpy import linalg

# Calcula a solução do sistema não-linear aplicando o método de Newton

def F(x, V1, V2, V_treshold, E1, E2, alpha1, alpha2):
    y = np.zeros(2)

    voltage = ((V2-V_treshold)/(V1-V_treshold))
    E1_coef = (E1**2 - 1)/E1
    E2_coef = (E2**2 - 1)/E2

    y[0] = -(2 + 2*voltage)*x[0]**2 + ( E1_coef * ((E2**2 + 1)/((1-E2**2)*np.sqrt(alpha2))) + (2+voltage)*((E1**2 + 1)/E1) )*x[0] +
                                                                                  ( E1_coef*(-2*E2)/((1-E2**2)*np.sqrt(alpha2)) )*x[0]*x[1] - 2

    y[1] = ( (E2**2 + 1)/(E2*np.sqrt(alpha2)) - voltage*E2_coef*((E2**2 + 1)/(1-E1**2)) )*x[1] +
                                                                                  ( voltage/(1-E1**2)*E2_coef )*x[0]*x[1] - 2/np.sqrt(alpha2)

    return y

def JF(x, V1, V2, V_treshold, E1, E2, alpha1, alpha2):

    y = np.zeros((2,2))

    y[0,0] = -2*(2 + 2*voltage)*x[0] + ( E1_coef * ((E2**2 + 1)/((1-E2**2)*np.sqrt(alpha2))) + (2+voltage)*((E1**2 + 1)/E1) ) +
                                                                                  ( E1_coef*(-2*E2)/((1-E2**2)*np.sqrt(alpha2)) )*x[1]
    y[0,1] = ( E1_coef*(-2*E2)/((1-E2**2)*np.sqrt(alpha2)) )*x[0]
    y[1,0] = ( voltage/(1-E1**2)*E2_coef )*x[1]
    y[1,1] = ( (E2**2 + 1)/(E2*np.sqrt(alpha2)) - voltage*E2_coef*((E2**2 + 1)/(1-E1**2)) ) + ( voltage/(1-E1**2)*E2_coef )*x[0]

    return y

def main(V1, V2, V_treshold, E1, E2, alpha1, alpha2):
    x = np.ones(2)*0.0000001

    for i in range(100):
        delta = -np.linalg.inv(JF(x, V1, V2, V_treshold, E1, E2, alpha1, alpha2)).dot(F(x, V1, V2, V_treshold, E1, E2, alpha1, alpha2))
        x = x + delta

    return x[0], x[1]
