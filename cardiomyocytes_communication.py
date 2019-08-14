#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Calcula a probabilidade de os connexons apresetarem os estados: HH, HL e LH

Utiliza as equações sugeridas por Kilinc, 2013

Autor: Ittalo dos Santos Silva <ittalosantoss@gmail.com>
"""

import numpy as np
from scipy.integrate import odeint
from cardiomyopy import *

def combination(N): # Cria uma matriz com todas as combinações possíveis de três números cuja soma seja igual a N

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
            raise Exception("Numero de combinacao p invalido")

    return combinations

def differential_equation(y, t, b1, a1, b2, a2): # Resolve equações de diferenciais do primeiro grau

    dydt1 = b1*(1 - (y[0] + y[1])) - a1*y[0]
    dydt2 = b2*(1 - (y[0] + y[1])) - a2*y[1]
    """
    # Condição de parada
    if y[0]>=1 and dydt1>=0:
        dhdt1 = 0
    if y[1]>=1 and dydt2>=0:
        dhdt2 = 0
    """
    dydt = [dydt1,dydt2]

    return dydt

def main(): # Função principal

    Rm = 2 # Resistividade da membrana do cardiomiócito [Ohm*m2]
    L = 100e-6 # Comprimento do cardiomiócito [m]
    r = 11e-6 # Raio do cardiomiócito [m]

    S = 2*np.pi()*r**2 + 2*np.pi()*r*L # Área da superfície da membrana do cardiomiócito [m2]
    V1 = -87e-3 # Potencial de repouso da membrana do cardiomiócito [V]
    V2 = 28.8e-3 # Potencial de equilíbrio do Sódio (Na+) [V]
    V_treshold = -34.8e-3 # Tensão limiar de ativação para liberar fluxo de Na+ [V]

    alpha1 = 1
    alpha2 = 108
    Ri = 2.5 # Resistividade citoplasmática média [Ohm*m]
    Re = 1 # Resistividade extracelular média [Ohm*m]
    Vi = np.pi()*(r**2)*L # Volume do citoplasma do cardiomiócito (intracelular) [m3]
    Ve = Vi/5 # Volume médio extracelular do cardiomiócito (extamente ao redor da célula) [m3]

    gamma = Ri/Vi + Re/Ve
    E1 = np.exp(np.sqrt((alpha1*gamma*S*L**2)/Rm))
    E2 = np.exp(np.sqrt((alpha2*gamma*S*L**2)/Rm))

    # Obtenção do mu1 e mu2 - Sistema de equações não-lineares do segundo grau
    mu1, mu2 = utils.main(V1, V2, V_treshold, E1, E2, alpha1, alpha2)

    # Resistência Crítica do GJ - Equação 2 do Kilinc 2013
    R_gj_critical = L*np.sqrt(gamma*Rm/S)*((E2**2 - 2*mu2*E2 + 1)/((1 - E2**2)*np.sqrt(alpha2)) +
                                             (V2 - V_treshold)/(V1 - V_treshold)*(E1**2 - 2*mu1*E1 + 1)/((1 - E1**2)*np.sqrt(alpha1)))

    # Condutância Crítica do GJ
    G_gj_critical = 1/(R_gj_critical)

    # Constantes para a proteína connexin43 do GJ
    lamb = 0.69
    A_alpha = 0.04e3 # [1/V]
    A_beta = 0.07e3 # [1/V]
    V0 = 62e-3 # Tensão da junção para satisfazer a igualdade de lamb nas equações [V]

    G_HH = 73e-12 # Condutância aberto-aberto
    G_HL = 12e-12 # Condutância aberto-fechado
    G_LH = 12e-12 # Condutância fechado-aberto

    N = np.arange(10, 20, 1, dtype = np.int32) # Número total de canais GJ
    Pt = Pr = np.zeros(N.size) # Potenciais de membrana do transmissor e do receptor

    # Variáveis de tempo
    dt = 0.01e-3 # Passo de tempo [s]
    tc = 1 # Quantidade de ciclos

    # Inicializar a rede de cardiomiócitos
    cardNet = CardiomyocyteNetwork(G_gj_critical, lamb, A_alpha, A_beta, V0, G_HH, G_HL, G_LH)
    cardNet.initialize_parameters(Pt, Pr)
    cardNet.set_time(dt, tc)

    Pf = np.zeros((num_curves, num_curves), dtype=np.float) # Probabilidade de falha da propagação
    N_plot = np.zeros((num_curves, num_curves)) # Vetor auxiliar para o gráfico
    Vj_plot = np.zeros((num_curves, num_curves)) # Vetor auxiliar para o gráfico

    for i in range(len(N)):

        for j in range(len(Vj)):
            N_plot[j, i] = N[i]
            Vj_plot[i, j] = Vj[j]

            # Constantes de transição dos estados
            alpha1_vj = lamb*math.exp(-A_alpha*(Vj[j] - V0))
            alpha2_vj = lamb*math.exp(A_alpha*(Vj[j] + V0))
            beta1_vj = lamb*math.exp(A_beta*(Vj[j] - V0))
            beta2_vj = lamb*math.exp(-A_beta*(Vj[j] + V0))

            # *** Equações Diferenciais para a as probabilidades - Começo ***

            #t = np.linspace(0, N[i], N[i]) # Vetor de tempo para a solução
            y0 = [0, 0] # Condições iniciais
            y = odeint(differential_equation, y0, Vj, args=(beta1_vj, alpha1_vj, beta2_vj, alpha2_vj)) # Resolve a equação diferencial
            P_LH = y[:, 0]
            P_HL = y[:, 1]

            P_HH = 1 - (P_LH + P_HL) # A soma das probabilidades é igual a 1

            # *** Equações Diferenciais para a as probabilidades - Fim ***

            comb = combination(N[i]) # Retorna matriz com todas as combinações de três números cuja soma é igual a N

            for n_HH, n_LH, n_HL in zip(comb[:, 0], comb[:, 1], comb[:, 2]):

                G_gj = n_HH*G_HH + n_LH*G_LH + n_HL*G_HL # Calcula condutância total para cada combinação de connexons

                if G_gj < G_gj_critical: # Se a condutância da combinação for menor que a condutância crítica
                    Pf[i, j] += math.factorial(N[i])*((P_HH[j]**n_HH)*(P_LH[j]**n_LH)*(P_HL[j]**n_HL))/(math.factorial(n_HH)*math.factorial(n_LH)*math.factorial(n_HL)) # Realiza a distribuição multinomial para a condição de falha na propagação

    # Gerar o gráfico das probabilidades
    probability_grapich(Vj_plot, N_plot, Pf, N)

if __name__ == "__main__":
    main()
