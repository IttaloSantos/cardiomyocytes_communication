# -*- coding: utf-8 -*-

try:
    from PySide import QtWidgets
except:
    from PyQt5 import QtWidgets


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

    def cardiomyocyte_transmitter():

    def cardiomyocyte_receiver():
