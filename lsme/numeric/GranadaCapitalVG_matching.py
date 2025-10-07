## Symbols used:
# ydbar
# gdVG
# yd
# yl
# guVG
# ylbar
# MVG
# gqVG
# muH
# lam
# yubar
# yu

import sys
import os
import numpy as np

import lsme.numeric.matchingresult as matchingresult

class GranadaVGMatchingResult(matchingresult.GenericMatchingResult):
    def __init__(self, name='VG'):
        super().__init__(name)
        self.MVG = 1
        self.gqVG = np.ones((3, 3))
        self.gqVGbar = np.ones((3, 3))
        self.guVG = np.ones((3, 3))
        self.guVGbar = np.ones((3, 3))
        self.gdVG = np.ones((3, 3))
        self.gdVGbar = np.ones((3, 3))
        self.exotic_params = ['MVG', 'gqVG', 'gqVGbar', 'guVG', 'guVGbar', 'gdVG', 'gdVGbar']
        self.nonvanishing = ['alphaOddbar', 'alphaOuubar', 'alphaOdd', 'alphaOqd8', 'alphaOqd8bar', 'alphaOqq1', 'alphaOqq1bar', 'alphaOqq3', 'alphaOqq3bar', 'alphaOqu8', 'alphaOqu8bar', 'alphaOud8', 'alphaOud8bar', 'alphaOuu']

    def alphaO3G(self, ):
        return 0

    def alphaO3Gt(self, ):
        return 0

    def alphaO3W(self, ):
        return 0

    def alphaO3Wt(self, ):
        return 0

    def alphaOH(self, ):
        return 0

    def alphaOHB(self, ):
        return 0

    def alphaOHBox(self, ):
        return 0

    def alphaOHBt(self, ):
        return 0

    def alphaOHD(self, ):
        return 0

    def alphaOHG(self, ):
        return 0

    def alphaOHGt(self, ):
        return 0

    def alphaOHW(self, ):
        return 0

    def alphaOHWB(self, ):
        return 0

    def alphaOHWBt(self, ):
        return 0

    def alphaOHWt(self, ):
        return 0

    def alphaOddbar(self, mif1,mif2,mif3,mif4):
        return (-1/2 * (self.MVG)**(-2) * self.gdVG[mif1,mif4] * self.gdVG[mif3,mif2] + 1/6 * (self.MVG)**(-2) * self.gdVG[mif1,mif2] * self.gdVG[mif3,mif4] + -1/2 * (self.MVG)**(-2) * self.gdVG[mif2,mif3] * self.gdVG[mif4,mif1] + 1/6 * (self.MVG)**(-2) * self.gdVG[mif2,mif1] * self.gdVG[mif4,mif3])

    def alphaOeebar(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOHdbar(self, mif1,mif2):
        return 0

    def alphaOHebar(self, mif1,mif2):
        return 0

    def alphaOHl1bar(self, mif1,mif2):
        return 0

    def alphaOHl3bar(self, mif1,mif2):
        return 0

    def alphaOHq1bar(self, mif1,mif2):
        return 0

    def alphaOHq3bar(self, mif1,mif2):
        return 0

    def alphaOHubar(self, mif1,mif2):
        return 0

    def alphaOllbar(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOuubar(self, mif1,mif2,mif3,mif4):
        return (-1/2 * (self.MVG)**(-2) * self.guVG[mif1,mif4] * self.guVG[mif3,mif2] + 1/6 * (self.MVG)**(-2) * self.guVG[mif1,mif2] * self.guVG[mif3,mif4] + -1/2 * (self.MVG)**(-2) * self.guVG[mif2,mif3] * self.guVG[mif4,mif1] + 1/6 * (self.MVG)**(-2) * self.guVG[mif2,mif1] * self.guVG[mif4,mif3])

    def alphaWeinberg(self, mif1,mif2):
        return 0

    def alphaWeinbergbar(self, mif1,mif2):
        return 0

    def alphaOdB(self, mif1,mif2):
        return 0

    def alphaOdBbar(self, mif1,mif2):
        return 0

    def alphaOdd(self, mif1,mif2,mif3,mif4):
        return (-1/2 * (self.MVG)**(-2) * self.gdVG[mif1,mif4] * self.gdVG[mif3,mif2] + 1/6 * (self.MVG)**(-2) * self.gdVG[mif1,mif2] * self.gdVG[mif3,mif4] + -1/2 * (self.MVG)**(-2) * self.gdVG[mif2,mif3] * self.gdVG[mif4,mif1] + 1/6 * (self.MVG)**(-2) * self.gdVG[mif2,mif1] * self.gdVG[mif4,mif3])

    def alphaOdG(self, mif1,mif2):
        return 0

    def alphaOdGbar(self, mif1,mif2):
        return 0

    def alphaOdH(self, mif1,mif2):
        return 0

    def alphaOdHbar(self, mif1,mif2):
        return 0

    def alphaOdW(self, mif1,mif2):
        return 0

    def alphaOdWbar(self, mif1,mif2):
        return 0

    def alphaOeB(self, mif1,mif2):
        return 0

    def alphaOeBbar(self, mif1,mif2):
        return 0

    def alphaOed(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOedbar(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOee(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOeH(self, mif1,mif2):
        return 0

    def alphaOeHbar(self, mif1,mif2):
        return 0

    def alphaOeu(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOeubar(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOeW(self, mif1,mif2):
        return 0

    def alphaOeWbar(self, mif1,mif2):
        return 0

    def alphaOHd(self, mif1,mif2):
        return 0

    def alphaOHe(self, mif1,mif2):
        return 0

    def alphaOHl1(self, mif1,mif2):
        return 0

    def alphaOHl3(self, mif1,mif2):
        return 0

    def alphaOHq1(self, mif1,mif2):
        return 0

    def alphaOHq3(self, mif1,mif2):
        return 0

    def alphaOHu(self, mif1,mif2):
        return 0

    def alphaOHud(self, mif1,mif2):
        return 0

    def alphaOHudbar(self, mif1,mif2):
        return 0

    def alphaOld(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOldbar(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOle(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOlebar(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOledq(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOledqbar(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOlequ1(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOlequ1bar(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOlequ3(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOlequ3bar(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOll(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOlq1(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOlq1bar(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOlq3(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOlq3bar(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOlu(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOlubar(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOqd1(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOqd1bar(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOqd8(self, mif1,mif2,mif3,mif4):
        return -4 * (self.MVG)**(-2) * self.gdVG[mif3,mif4] * self.gqVG[mif1,mif2]

    def alphaOqd8bar(self, mif1,mif2,mif3,mif4):
        return -4 * (self.MVG)**(-2) * self.gdVG[mif3,mif4] * self.gqVG[mif1,mif2]

    def alphaOqe(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOqebar(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOqq1(self, mif1,mif2,mif3,mif4):
        return (-1/4 * (self.MVG)**(-2) * self.gqVG[mif1,mif4] * self.gqVG[mif3,mif2] + 1/6 * (self.MVG)**(-2) * self.gqVG[mif1,mif2] * self.gqVG[mif3,mif4] + -1/4 * (self.MVG)**(-2) * self.gqVG[mif2,mif3] * self.gqVG[mif4,mif1] + 1/6 * (self.MVG)**(-2) * self.gqVG[mif2,mif1] * self.gqVG[mif4,mif3])

    def alphaOqq1bar(self, mif1,mif2,mif3,mif4):
        return (-1/4 * (self.MVG)**(-2) * self.gqVG[mif1,mif4] * self.gqVG[mif3,mif2] + 1/6 * (self.MVG)**(-2) * self.gqVG[mif1,mif2] * self.gqVG[mif3,mif4] + -1/4 * (self.MVG)**(-2) * self.gqVG[mif2,mif3] * self.gqVG[mif4,mif1] + 1/6 * (self.MVG)**(-2) * self.gqVG[mif2,mif1] * self.gqVG[mif4,mif3])

    def alphaOqq3(self, mif1,mif2,mif3,mif4):
        return (-1/4 * (self.MVG)**(-2) * self.gqVG[mif1,mif4] * self.gqVG[mif3,mif2] + -1/4 * (self.MVG)**(-2) * self.gqVG[mif2,mif3] * self.gqVG[mif4,mif1])

    def alphaOqq3bar(self, mif1,mif2,mif3,mif4):
        return (-1/4 * (self.MVG)**(-2) * self.gqVG[mif1,mif4] * self.gqVG[mif3,mif2] + -1/4 * (self.MVG)**(-2) * self.gqVG[mif2,mif3] * self.gqVG[mif4,mif1])

    def alphaOqu1(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOqu1bar(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOqu8(self, mif1,mif2,mif3,mif4):
        return -4 * (self.MVG)**(-2) * self.gqVG[mif1,mif2] * self.guVG[mif3,mif4]

    def alphaOqu8bar(self, mif1,mif2,mif3,mif4):
        return -4 * (self.MVG)**(-2) * self.gqVG[mif1,mif2] * self.guVG[mif3,mif4]

    def alphaOquqd1(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOquqd1bar(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOquqd8(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOquqd8bar(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOuB(self, mif1,mif2):
        return 0

    def alphaOuBbar(self, mif1,mif2):
        return 0

    def alphaOud1(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOud1bar(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOud8(self, mif1,mif2,mif3,mif4):
        return -4 * (self.MVG)**(-2) * self.gdVG[mif3,mif4] * self.guVG[mif1,mif2]

    def alphaOud8bar(self, mif1,mif2,mif3,mif4):
        return -4 * (self.MVG)**(-2) * self.gdVG[mif3,mif4] * self.guVG[mif1,mif2]

    def alphaOuG(self, mif1,mif2):
        return 0

    def alphaOuGbar(self, mif1,mif2):
        return 0

    def alphaOuH(self, mif1,mif2):
        return 0

    def alphaOuHbar(self, mif1,mif2):
        return 0

    def alphaOuu(self, mif1,mif2,mif3,mif4):
        return (-1/2 * (self.MVG)**(-2) * self.guVG[mif1,mif4] * self.guVG[mif3,mif2] + 1/6 * (self.MVG)**(-2) * self.guVG[mif1,mif2] * self.guVG[mif3,mif4] + -1/2 * (self.MVG)**(-2) * self.guVG[mif2,mif3] * self.guVG[mif4,mif1] + 1/6 * (self.MVG)**(-2) * self.guVG[mif2,mif1] * self.guVG[mif4,mif3])

    def alphaOuW(self, mif1,mif2):
        return 0

    def alphaOuWbar(self, mif1,mif2):
        return 0
