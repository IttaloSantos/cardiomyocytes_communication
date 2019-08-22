#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Calcula a probabilidade de os connexons apresetarem os estados: HH, HL e LH

Utiliza as equações sugeridas por Kilinc, 2013

Autor: Ittalo dos Santos Silva <ittalosantoss@gmail.com>
"""

import numpy as np
from cardiomyopy import *

def main(): # Função principal

    Rm = 2 # Resistividade da membrana do cardiomiócito [Ohm*m2]
    L = 100e-6 # Comprimento do cardiomiócito [m]
    r = 11e-6 # Raio do cardiomiócito [m]

    S = 2*np.pi*r**2 + 2*np.pi*r*L # Área da superfície da membrana do cardiomiócito [m2]
    V1 = -87#e-3 # Potencial de repouso da membrana do cardiomiócito [mV]
    V2 = 28.8#e-3 # Potencial de equilíbrio do Sódio (Na+) [mV]
    V_treshold = -34.8#e-3 # Tensão limiar de ativação para liberar fluxo de Na+ [mV]

    alpha1 = 1
    alpha2 = 108
    Ri = 2.5 # Resistividade citoplasmática média [Ohm*m]
    Re = 1 # Resistividade extracelular média [Ohm*m]
    Vi = np.pi*(r**2)*L # Volume do citoplasma do cardiomiócito (intracelular) [m3]
    Ve = Vi/5 # Volume médio extracelular do cardiomiócito (extamente ao redor da célula) [m3]

    gamma = Ri/Vi + Re/Ve
    E1 = np.exp(np.sqrt((alpha1*gamma*S*L**2)/Rm))
    E2 = np.exp(np.sqrt((alpha2*gamma*S*L**2)/Rm))

    # Obtenção do mu1 e mu2 - Sistema de equações não-lineares do segundo grau
    mu1, mu2 = utils.non_linear_equations_solve(V1, V2, V_treshold, E1, E2, alpha1, alpha2)

    # Resistência Crítica do GJ - Equação 2 do Kilinc 2013
    R_gj_critical = L*np.sqrt(gamma*Rm/S)*((E2**2 - 2*mu2*E2 + 1)/((1 - E2**2)*np.sqrt(alpha2)) + (V2 - V_treshold)/(V1 - V_treshold)*(E1**2 - 2*mu1*E1 + 1)/((1 - E1**2)*np.sqrt(alpha1)))

    # Condutância Crítica do GJ
    G_gj_critical = 1/(R_gj_critical)

    # Constantes para a proteína connexin43 do GJ
    lamb = 0.69
    A_alpha = 0.04 # [1/mV]
    A_beta = 0.07 # [1/mV]
    V0 = 62e-3 # Tensão da junção para satisfazer a igualdade de lamb nas equações [V]

    G_HH = 73e-12 # Condutância aberto-aberto
    G_HL = 12e-12 # Condutância aberto-fechado
    G_LH = 12e-12 # Condutância fechado-aberto

    N = np.arange(100, 200, 1, dtype = np.int32) # Número total de canais GJ
    Pt = Pr = np.zeros(N.size) # Potenciais de membrana do transmissor e do receptor
    Pt_max = 40 # máximo valor do potencial de membrana do trasmissor
    Pr_max = 0 # máximo valor do potencial de membrana do receptor
    Vj = Pt-Pr # Tensão da junção - Diferença entre os potenciais de membrana do Tc e do Rc

    # Variáveis de tempo
    dt = 0.01e-3 # Passo de tempo [s]
    tc = 1 # Quantidade de ciclos

    # Inicializar a rede de cardiomiócitos
    cardNet = CardiomyocyteNetwork(G_gj_critical, lamb, A_alpha, A_beta, V0, G_HH, G_HL, G_LH)
    cardNet.initialize_parameters(Pt, Pr, Vj)
    cardNet.set_time(dt, tc)

    # Iniciar a comunicação
    cardNet.communication(N, Pt_max, Pr_max)

    # Gerar o gráfico das probabilidades
    cardNet.probability_graphic()

if __name__ == "__main__":
    main()
