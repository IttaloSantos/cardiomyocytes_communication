# -*- coding: utf-8 -*-

import numpy as np
import math as math
from time import clock
from decimal import Decimal

from cardiomyopy import utils as utils
from cardiomyopy.cardiomyocyte import Cardiomyocyte

from os import makedirs
from os.path import exists

class CardiomyocyteNetwork(object):
    """
    Classe representando uma rede de cardiomiócitos

    param G_gj_critical: condutância crítica
    param lamb: lambda
    param A_alpha: constante para o cálculo da condutância instantânea [1/V]
    param A_beta: constante para o cálculo da condutância instantânea [1/V]
    param V0: Tensão da junção para satisfazer a igualdade de lamb nas equações [V]
    param G_HH, G_HL, G_LH: condutância dos respectivos estados [Siemens]
    """

    def __init__(self, G_gj_critical, lamb, A_alpha, A_beta, V0, G_HH, G_HL, G_LH):
        # Construtor da rede

        self.G_gj_critical = G_gj_critical
        self.lamb = lamb
        self.A_alpha = A_alpha
        self.A_beta = A_beta
        self.V0 = V0
        self.G_HH = G_HH
        self.G_HL = G_HL
        self.G_LH = G_LH

    def initialize_parameters(self, Pt, Pr, Vj):
        # param Pt, Pr: potenciais de membrana do transmissor e do receptor
        # param Vj: Tensão de junção - Diferença entre os potenciais de membrana do Tc e do Rc

        self.Vj = Vj
        self.cardiomyocytes = Cardiomyocyte(Pt, Pr)

    def set_time(self, dt, tc):
        # param dt: passo de tempo
        # param tc: quantidade de ciclos

        self.T = 0
        self.dt = dt
        self.t0 = clock()
        self.tc = tc

    def time_step(self):
        # Incremento no tempo do processo
        # param T: variável de tempo

        self.T = clock() - self.t0

    def set_Vj(self):
        # Retorna a diferença entre os potenciais de membrana do TC e RC
        self.Vj = self.cardiomyocytes.DDP()

    def probability_graphic(self):
        # Salva os resultados em um arquivo CSV e gera o gráfico
        utils.dump_results(self.Pf, self.N, self.Vj)

    def failure_probability_calc(self, i, j, N, n_HH, n_LH, n_HL):
        # Realiza a distribuição multinomial para a condição de falha na propagação

        self.Pf[i, j] += float(math.factorial(N)*Decimal((self.P_HH**n_HH)*(self.P_LH**n_LH)*(self.P_HL**n_HL))/(math.factorial(n_HH)*math.factorial(n_LH)*math.factorial(n_HL)))

    def communication(self, N, Pt_max, Pr_max):
        # Estabelece a comunicação propriamente dita

        # param Pt_max, Pr_max: máximos valores dos potenciais de membrana
        # param N: número total de canais GJ

        self.cardiomyocytes.stimulate(Pt_max, Pr_max, N.size)
        self.set_Vj()
        self.N = N

        self.Pf = np.zeros((N.size, N.size), dtype=np.float) # Probabilidade de falha da propagação
        for i in range(N.size):
            self.time_step()
            print("Tempo do processo: ", round(self.T, 1), "s")
            for j in range(self.Vj.size):
                # Constantes de transição dos estados
                alpha1_vj = self.lamb*np.exp(-self.A_alpha*(self.Vj[j] - self.V0))
                alpha2_vj = self.lamb*np.exp(self.A_alpha*(self.Vj[j] + self.V0))
                beta1_vj = self.lamb*np.exp(self.A_beta*(self.Vj[j] - self.V0))
                beta2_vj = self.lamb*np.exp(-self.A_beta*(self.Vj[j] + self.V0))

                # *** Equações Diferenciais para as probabilidades - Começo ***
                t = np.arange(0, 1500.01e-3, 0.01e-3)
                p0 = [0, 0] # Variável de probabilidade: Condições iniciais
                p = utils.differential_equation_solve(p0, t, beta1_vj, alpha1_vj, beta2_vj, alpha2_vj) # Variável de probabilidade: solução das equações diferenciais

                self.P_LH = np.mean(p[:, 0]) # Probabilidade do estado Low-High dos canais GJ - Média temporal
                self.P_HL = np.mean(p[:, 1]) # Probabilidade do estado High-Low dos canais GJ - Média temporal
                self.P_HH = 1 - (self.P_LH + self.P_HL) # Probabilidade do estado High-High dos canais GJ, ignorando o estado Low-Low

                # *** Equações Diferenciais para a as probabilidades - Fim ***
                comb = utils.combination(N[i]) # Retorna matriz com todas as combinações de três números cuja soma é igual a N

                for n_HH, n_LH, n_HL in zip(comb[:, 0], comb[:, 1], comb[:, 2]):

                    G_gj = n_HH*self.G_HH + n_LH*self.G_LH + n_HL*self.G_HL # Calcula condutância total para cada combinação de connexons

                    if G_gj > self.G_gj_critical: # Se a condutância da combinação for menor que a condutância crítica
                        self.failure_probability_calc(i, j, N[i], n_HH, n_LH, n_HL)
