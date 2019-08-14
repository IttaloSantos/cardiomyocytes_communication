# -*- coding: utf-8 -*-

import numpy as np
import math as math

from cardiomyopy import utils as utils
from cardiomyopy.cardiomyocyte import Cardiomyocyte

from os import makedirs
from os.path import exists

try:
    from PySide import QtWidgets
except:
    from PyQt5 import QtWidgets


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
        self.G_HH = G_HH
        self.G_HL = G_HL
        self.G_LH = G_LH

    def initialize_parameters(Pt, Pr):
        # param Vj: Tensão da junção - Valor de tensão do vetor Vj do Kilinc 2013
        # param Pt, Pr: potenciais de membrana do transmissor e do receptor

        self.cardiomyocytes = Cardiomyocytes(Pt, Pr)
