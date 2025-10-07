## Symbols used:
# ydbar
# yd
# glVW
# yl
# ylbar
# muH
# MVW
# lam
# gphiVWbar
# gphiVW
# yubar
# gqVW
# yu

import sys
import os
import numpy as np

import lsme.numeric.matchingresult as matchingresult

class GranadaVWMatchingResult(matchingresult.GenericMatchingResult):
    def __init__(self, name='VW'):
        super().__init__(name)
        self.MVW = 1
        self.glVW = np.ones((3, 3))
        self.glVWbar = np.ones((3, 3))
        self.gqVW = np.ones((3, 3))
        self.gqVWbar = np.ones((3, 3))
        self.gphiVW = 1
        self.gphiVWbar = 1
        self.exotic_params = ['MVW', 'glVW', 'glVWbar', 'gqVW', 'gqVWbar', 'gphiVW', 'gphiVWbar']
        self.nonvanishing = ['alphaOH', 'alphaOHBox', 'alphaOHD', 'alphaOHl3bar', 'alphaOHq3bar', 'alphaOllbar', 'alphaOdH', 'alphaOdHbar', 'alphaOeH', 'alphaOeHbar', 'alphaOHl3', 'alphaOHq3', 'alphaOll', 'alphaOlq3', 'alphaOlq3bar', 'alphaOqq3', 'alphaOqq3bar', 'alphaOuH', 'alphaOuHbar']

    def alphaO3G(self, ):
        return 0

    def alphaO3Gt(self, ):
        return 0

    def alphaO3W(self, ):
        return 0

    def alphaO3Wt(self, ):
        return 0

    def alphaOH(self, ):
        return -1 * self.gphiVW * self.gphiVWbar * self.lam * (self.MVW)**(-2)

    def alphaOHB(self, ):
        return 0

    def alphaOHBox(self, ):
        return (-1/16 * (self.gphiVW)**(2) * (self.MVW)**(-2) + -1/4 * self.gphiVW * self.gphiVWbar * (self.MVW)**(-2) + -1/16 * (self.gphiVWbar)**(2) * (self.MVW)**(-2))

    def alphaOHBt(self, ):
        return 0

    def alphaOHD(self, ):
        return (-1/8 * (self.gphiVW)**(2) * (self.MVW)**(-2) + 1/4 * self.gphiVW * self.gphiVWbar * (self.MVW)**(-2) + -1/8 * (self.gphiVWbar)**(2) * (self.MVW)**(-2))

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
        return 0

    def alphaOeebar(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOHdbar(self, mif1,mif2):
        return 0

    def alphaOHebar(self, mif1,mif2):
        return 0

    def alphaOHl1bar(self, mif1,mif2):
        return 0

    def alphaOHl3bar(self, mif1,mif2):
        return (-1/8 * self.gphiVW * (self.MVW)**(-2) * self.glVW[mif1,mif2] + -1/8 * self.gphiVWbar * (self.MVW)**(-2) * self.glVW[mif1,mif2])

    def alphaOHq1bar(self, mif1,mif2):
        return 0

    def alphaOHq3bar(self, mif1,mif2):
        return (-1/8 * self.gphiVW * (self.MVW)**(-2) * self.gqVW[mif1,mif2] + -1/8 * self.gphiVWbar * (self.MVW)**(-2) * self.gqVW[mif1,mif2])

    def alphaOHubar(self, mif1,mif2):
        return 0

    def alphaOllbar(self, mif1,mif2,mif3,mif4):
        return (-1/8 * (self.MVW)**(-2) * self.glVW[mif1,mif4] * self.glVW[mif3,mif2] + 1/16 * (self.MVW)**(-2) * self.glVW[mif1,mif2] * self.glVW[mif3,mif4] + -1/8 * (self.MVW)**(-2) * self.glVW[mif2,mif3] * self.glVW[mif4,mif1] + 1/16 * (self.MVW)**(-2) * self.glVW[mif2,mif1] * self.glVW[mif4,mif3])

    def alphaOuubar(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaWeinberg(self, mif1,mif2):
        return 0

    def alphaWeinbergbar(self, mif1,mif2):
        return 0

    def alphaOdB(self, mif1,mif2):
        return 0

    def alphaOdBbar(self, mif1,mif2):
        return 0

    def alphaOdd(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOdG(self, mif1,mif2):
        return 0

    def alphaOdGbar(self, mif1,mif2):
        return 0

    def alphaOdH(self, mif1,mif2):
        return (1/16 * self.gphiVW * (self.MVW)**(-2) * self.gqVW[0,mif1] * self.yd[0,mif2] + -1/16 * self.gphiVWbar * (self.MVW)**(-2) * self.gqVW[0,mif1] * self.yd[0,mif2] + 1/16 * self.gphiVW * (self.MVW)**(-2) * self.gqVW[mif1,0] * self.yd[0,mif2] + -1/16 * self.gphiVWbar * (self.MVW)**(-2) * self.gqVW[mif1,0] * self.yd[0,mif2] + 1/16 * self.gphiVW * (self.MVW)**(-2) * self.gqVW[1,mif1] * self.yd[1,mif2] + -1/16 * self.gphiVWbar * (self.MVW)**(-2) * self.gqVW[1,mif1] * self.yd[1,mif2] + 1/16 * self.gphiVW * (self.MVW)**(-2) * self.gqVW[mif1,1] * self.yd[1,mif2] + -1/16 * self.gphiVWbar * (self.MVW)**(-2) * self.gqVW[mif1,1] * self.yd[1,mif2] + 1/16 * self.gphiVW * (self.MVW)**(-2) * self.gqVW[2,mif1] * self.yd[2,mif2] + -1/16 * self.gphiVWbar * (self.MVW)**(-2) * self.gqVW[2,mif1] * self.yd[2,mif2] + 1/16 * self.gphiVW * (self.MVW)**(-2) * self.gqVW[mif1,2] * self.yd[2,mif2] + -1/16 * self.gphiVWbar * (self.MVW)**(-2) * self.gqVW[mif1,2] * self.yd[2,mif2] + -1/16 * (self.gphiVW)**(2) * (self.MVW)**(-2) * self.yd[mif1,mif2] + -1/4 * self.gphiVW * self.gphiVWbar * (self.MVW)**(-2) * self.yd[mif1,mif2] + 1/16 * (self.gphiVWbar)**(2) * (self.MVW)**(-2) * self.yd[mif1,mif2])

    def alphaOdHbar(self, mif1,mif2):
        return (-1/16 * self.gphiVW * (self.MVW)**(-2) * self.gqVW[0,mif1] * self.ydbar[0,mif2] + 1/16 * self.gphiVWbar * (self.MVW)**(-2) * self.gqVW[0,mif1] * self.ydbar[0,mif2] + -1/16 * self.gphiVW * (self.MVW)**(-2) * self.gqVW[mif1,0] * self.ydbar[0,mif2] + 1/16 * self.gphiVWbar * (self.MVW)**(-2) * self.gqVW[mif1,0] * self.ydbar[0,mif2] + -1/16 * self.gphiVW * (self.MVW)**(-2) * self.gqVW[1,mif1] * self.ydbar[1,mif2] + 1/16 * self.gphiVWbar * (self.MVW)**(-2) * self.gqVW[1,mif1] * self.ydbar[1,mif2] + -1/16 * self.gphiVW * (self.MVW)**(-2) * self.gqVW[mif1,1] * self.ydbar[1,mif2] + 1/16 * self.gphiVWbar * (self.MVW)**(-2) * self.gqVW[mif1,1] * self.ydbar[1,mif2] + -1/16 * self.gphiVW * (self.MVW)**(-2) * self.gqVW[2,mif1] * self.ydbar[2,mif2] + 1/16 * self.gphiVWbar * (self.MVW)**(-2) * self.gqVW[2,mif1] * self.ydbar[2,mif2] + -1/16 * self.gphiVW * (self.MVW)**(-2) * self.gqVW[mif1,2] * self.ydbar[2,mif2] + 1/16 * self.gphiVWbar * (self.MVW)**(-2) * self.gqVW[mif1,2] * self.ydbar[2,mif2] + 1/16 * (self.gphiVW)**(2) * (self.MVW)**(-2) * self.ydbar[mif1,mif2] + -1/4 * self.gphiVW * self.gphiVWbar * (self.MVW)**(-2) * self.ydbar[mif1,mif2] + -1/16 * (self.gphiVWbar)**(2) * (self.MVW)**(-2) * self.ydbar[mif1,mif2])

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
        return (1/16 * self.gphiVW * (self.MVW)**(-2) * self.glVW[0,mif1] * self.yl[0,mif2] + -1/16 * self.gphiVWbar * (self.MVW)**(-2) * self.glVW[0,mif1] * self.yl[0,mif2] + 1/16 * self.gphiVW * (self.MVW)**(-2) * self.glVW[mif1,0] * self.yl[0,mif2] + -1/16 * self.gphiVWbar * (self.MVW)**(-2) * self.glVW[mif1,0] * self.yl[0,mif2] + 1/16 * self.gphiVW * (self.MVW)**(-2) * self.glVW[1,mif1] * self.yl[1,mif2] + -1/16 * self.gphiVWbar * (self.MVW)**(-2) * self.glVW[1,mif1] * self.yl[1,mif2] + 1/16 * self.gphiVW * (self.MVW)**(-2) * self.glVW[mif1,1] * self.yl[1,mif2] + -1/16 * self.gphiVWbar * (self.MVW)**(-2) * self.glVW[mif1,1] * self.yl[1,mif2] + 1/16 * self.gphiVW * (self.MVW)**(-2) * self.glVW[2,mif1] * self.yl[2,mif2] + -1/16 * self.gphiVWbar * (self.MVW)**(-2) * self.glVW[2,mif1] * self.yl[2,mif2] + 1/16 * self.gphiVW * (self.MVW)**(-2) * self.glVW[mif1,2] * self.yl[2,mif2] + -1/16 * self.gphiVWbar * (self.MVW)**(-2) * self.glVW[mif1,2] * self.yl[2,mif2] + -1/16 * (self.gphiVW)**(2) * (self.MVW)**(-2) * self.yl[mif1,mif2] + -1/4 * self.gphiVW * self.gphiVWbar * (self.MVW)**(-2) * self.yl[mif1,mif2] + 1/16 * (self.gphiVWbar)**(2) * (self.MVW)**(-2) * self.yl[mif1,mif2])

    def alphaOeHbar(self, mif1,mif2):
        return (-1/16 * self.gphiVW * (self.MVW)**(-2) * self.glVW[0,mif1] * self.ylbar[0,mif2] + 1/16 * self.gphiVWbar * (self.MVW)**(-2) * self.glVW[0,mif1] * self.ylbar[0,mif2] + -1/16 * self.gphiVW * (self.MVW)**(-2) * self.glVW[mif1,0] * self.ylbar[0,mif2] + 1/16 * self.gphiVWbar * (self.MVW)**(-2) * self.glVW[mif1,0] * self.ylbar[0,mif2] + -1/16 * self.gphiVW * (self.MVW)**(-2) * self.glVW[1,mif1] * self.ylbar[1,mif2] + 1/16 * self.gphiVWbar * (self.MVW)**(-2) * self.glVW[1,mif1] * self.ylbar[1,mif2] + -1/16 * self.gphiVW * (self.MVW)**(-2) * self.glVW[mif1,1] * self.ylbar[1,mif2] + 1/16 * self.gphiVWbar * (self.MVW)**(-2) * self.glVW[mif1,1] * self.ylbar[1,mif2] + -1/16 * self.gphiVW * (self.MVW)**(-2) * self.glVW[2,mif1] * self.ylbar[2,mif2] + 1/16 * self.gphiVWbar * (self.MVW)**(-2) * self.glVW[2,mif1] * self.ylbar[2,mif2] + -1/16 * self.gphiVW * (self.MVW)**(-2) * self.glVW[mif1,2] * self.ylbar[2,mif2] + 1/16 * self.gphiVWbar * (self.MVW)**(-2) * self.glVW[mif1,2] * self.ylbar[2,mif2] + 1/16 * (self.gphiVW)**(2) * (self.MVW)**(-2) * self.ylbar[mif1,mif2] + -1/4 * self.gphiVW * self.gphiVWbar * (self.MVW)**(-2) * self.ylbar[mif1,mif2] + -1/16 * (self.gphiVWbar)**(2) * (self.MVW)**(-2) * self.ylbar[mif1,mif2])

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
        return (-1/8 * self.gphiVW * (self.MVW)**(-2) * self.glVW[mif1,mif2] + -1/8 * self.gphiVWbar * (self.MVW)**(-2) * self.glVW[mif1,mif2])

    def alphaOHq1(self, mif1,mif2):
        return 0

    def alphaOHq3(self, mif1,mif2):
        return (-1/8 * self.gphiVW * (self.MVW)**(-2) * self.gqVW[mif1,mif2] + -1/8 * self.gphiVWbar * (self.MVW)**(-2) * self.gqVW[mif1,mif2])

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
        return (-1/8 * (self.MVW)**(-2) * self.glVW[mif1,mif4] * self.glVW[mif3,mif2] + 1/16 * (self.MVW)**(-2) * self.glVW[mif1,mif2] * self.glVW[mif3,mif4] + -1/8 * (self.MVW)**(-2) * self.glVW[mif2,mif3] * self.glVW[mif4,mif1] + 1/16 * (self.MVW)**(-2) * self.glVW[mif2,mif1] * self.glVW[mif4,mif3])

    def alphaOlq1(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOlq1bar(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOlq3(self, mif1,mif2,mif3,mif4):
        return -1/4 * (self.MVW)**(-2) * self.glVW[mif1,mif2] * self.gqVW[mif3,mif4]

    def alphaOlq3bar(self, mif1,mif2,mif3,mif4):
        return -1/4 * (self.MVW)**(-2) * self.glVW[mif1,mif2] * self.gqVW[mif3,mif4]

    def alphaOlu(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOlubar(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOqd1(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOqd1bar(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOqd8(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOqd8bar(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOqe(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOqebar(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOqq1(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOqq1bar(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOqq3(self, mif1,mif2,mif3,mif4):
        return -1/8 * (self.MVW)**(-2) * self.gqVW[mif1,mif2] * self.gqVW[mif3,mif4]

    def alphaOqq3bar(self, mif1,mif2,mif3,mif4):
        return -1/8 * (self.MVW)**(-2) * self.gqVW[mif1,mif2] * self.gqVW[mif3,mif4]

    def alphaOqu1(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOqu1bar(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOqu8(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOqu8bar(self, mif1,mif2,mif3,mif4):
        return 0

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
        return 0

    def alphaOud8bar(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOuG(self, mif1,mif2):
        return 0

    def alphaOuGbar(self, mif1,mif2):
        return 0

    def alphaOuH(self, mif1,mif2):
        return (-1/16 * self.gphiVW * (self.MVW)**(-2) * self.gqVW[0,mif1] * self.yu[0,mif2] + 1/16 * self.gphiVWbar * (self.MVW)**(-2) * self.gqVW[0,mif1] * self.yu[0,mif2] + -1/16 * self.gphiVW * (self.MVW)**(-2) * self.gqVW[mif1,0] * self.yu[0,mif2] + 1/16 * self.gphiVWbar * (self.MVW)**(-2) * self.gqVW[mif1,0] * self.yu[0,mif2] + -1/16 * self.gphiVW * (self.MVW)**(-2) * self.gqVW[1,mif1] * self.yu[1,mif2] + 1/16 * self.gphiVWbar * (self.MVW)**(-2) * self.gqVW[1,mif1] * self.yu[1,mif2] + -1/16 * self.gphiVW * (self.MVW)**(-2) * self.gqVW[mif1,1] * self.yu[1,mif2] + 1/16 * self.gphiVWbar * (self.MVW)**(-2) * self.gqVW[mif1,1] * self.yu[1,mif2] + -1/16 * self.gphiVW * (self.MVW)**(-2) * self.gqVW[2,mif1] * self.yu[2,mif2] + 1/16 * self.gphiVWbar * (self.MVW)**(-2) * self.gqVW[2,mif1] * self.yu[2,mif2] + -1/16 * self.gphiVW * (self.MVW)**(-2) * self.gqVW[mif1,2] * self.yu[2,mif2] + 1/16 * self.gphiVWbar * (self.MVW)**(-2) * self.gqVW[mif1,2] * self.yu[2,mif2] + 1/16 * (self.gphiVW)**(2) * (self.MVW)**(-2) * self.yu[mif1,mif2] + -1/4 * self.gphiVW * self.gphiVWbar * (self.MVW)**(-2) * self.yu[mif1,mif2] + -1/16 * (self.gphiVWbar)**(2) * (self.MVW)**(-2) * self.yu[mif1,mif2])

    def alphaOuHbar(self, mif1,mif2):
        return (1/16 * self.gphiVW * (self.MVW)**(-2) * self.gqVW[0,mif1] * self.yubar[0,mif2] + -1/16 * self.gphiVWbar * (self.MVW)**(-2) * self.gqVW[0,mif1] * self.yubar[0,mif2] + 1/16 * self.gphiVW * (self.MVW)**(-2) * self.gqVW[mif1,0] * self.yubar[0,mif2] + -1/16 * self.gphiVWbar * (self.MVW)**(-2) * self.gqVW[mif1,0] * self.yubar[0,mif2] + 1/16 * self.gphiVW * (self.MVW)**(-2) * self.gqVW[1,mif1] * self.yubar[1,mif2] + -1/16 * self.gphiVWbar * (self.MVW)**(-2) * self.gqVW[1,mif1] * self.yubar[1,mif2] + 1/16 * self.gphiVW * (self.MVW)**(-2) * self.gqVW[mif1,1] * self.yubar[1,mif2] + -1/16 * self.gphiVWbar * (self.MVW)**(-2) * self.gqVW[mif1,1] * self.yubar[1,mif2] + 1/16 * self.gphiVW * (self.MVW)**(-2) * self.gqVW[2,mif1] * self.yubar[2,mif2] + -1/16 * self.gphiVWbar * (self.MVW)**(-2) * self.gqVW[2,mif1] * self.yubar[2,mif2] + 1/16 * self.gphiVW * (self.MVW)**(-2) * self.gqVW[mif1,2] * self.yubar[2,mif2] + -1/16 * self.gphiVWbar * (self.MVW)**(-2) * self.gqVW[mif1,2] * self.yubar[2,mif2] + -1/16 * (self.gphiVW)**(2) * (self.MVW)**(-2) * self.yubar[mif1,mif2] + -1/4 * self.gphiVW * self.gphiVWbar * (self.MVW)**(-2) * self.yubar[mif1,mif2] + 1/16 * (self.gphiVWbar)**(2) * (self.MVW)**(-2) * self.yubar[mif1,mif2])

    def alphaOuu(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOuW(self, mif1,mif2):
        return 0

    def alphaOuWbar(self, mif1,mif2):
        return 0
