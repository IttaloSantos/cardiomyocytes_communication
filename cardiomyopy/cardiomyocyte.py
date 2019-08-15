# -*- coding: utf-8 -*-

import numpy as np

class Cardiomyocyte(object):
    """
    Classe que representa o cardiomiócito transmissor e o receptor.

    param Pt: potencial de membrana do cardiomiócito transmissor
    param Pr: potencial de membrana do cardiomiócito receptor
    """

    def __init__(self, Pt, Pr):
        # Construtor do cardiomiócito

        self.Pt = Pt
        self.Pr = Pr

    def stimulate(self, Pt_max, Pr_max, Nsize):
        # Dá o estímulo inicial aos cardiomiócitos para começar a comunicação

        # param Pt_max, Pr_max: máximos valores dos potenciais de membrana
        # param Nsize: tamanho do vetor número de canais GJ

        self.Pt = np.linspace(0, Pt_max, Nsize)
        if Pr_max != 0:
            self.Pr = np.linspace(0, Pr_max, Nsize)

    def DDP(self):
        # Retorna a diferença entre os potenciais de membrana do TC e RC

        return (self.Pt - self.Pr)
    """
    def cardiomyocyte_transmitter():

    def cardiomyocyte_receiver():
    """
