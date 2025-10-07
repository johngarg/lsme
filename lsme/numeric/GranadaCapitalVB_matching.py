## Symbols used:
# ydbar
# yd
# gphiVBbar
# gphiVB
# MVB
# yl
# ylbar
# muH
# guVB
# geVB
# gdVB
# gqVB
# lam
# yubar
# glVB
# yu

import sys
import os
import numpy as np

import lsme.numeric.matchingresult as matchingresult

class GranadaVBMatchingResult(matchingresult.GenericMatchingResult):
    def __init__(self, name='VB'):
        super().__init__(name)
        self.MVB = 1
        self.glVB = np.ones((3, 3))
        self.glVBbar = np.ones((3, 3))
        self.gqVB = np.ones((3, 3))
        self.gqVBbar = np.ones((3, 3))
        self.geVB = np.ones((3, 3))
        self.geVBbar = np.ones((3, 3))
        self.guVB = np.ones((3, 3))
        self.guVBbar = np.ones((3, 3))
        self.gdVB = np.ones((3, 3))
        self.gdVBbar = np.ones((3, 3))
        self.gphiVB = 1
        self.gphiVBbar = 1
        self.exotic_params = ['MVB', 'glVB', 'glVBbar', 'gqVB', 'gqVBbar', 'geVB', 'geVBbar', 'guVB', 'guVBbar', 'gdVB', 'gdVBbar', 'gphiVB', 'gphiVBbar']
        self.nonvanishing = ['alphaOHBox', 'alphaOHD', 'alphaOddbar', 'alphaOeebar', 'alphaOHdbar', 'alphaOHebar', 'alphaOHl1bar', 'alphaOHq1bar', 'alphaOHubar', 'alphaOllbar', 'alphaOuubar', 'alphaOdd', 'alphaOdH', 'alphaOdHbar', 'alphaOed', 'alphaOedbar', 'alphaOee', 'alphaOeH', 'alphaOeHbar', 'alphaOeu', 'alphaOeubar', 'alphaOHd', 'alphaOHe', 'alphaOHl1', 'alphaOHq1', 'alphaOHu', 'alphaOld', 'alphaOldbar', 'alphaOle', 'alphaOlebar', 'alphaOll', 'alphaOlq1', 'alphaOlq1bar', 'alphaOlu', 'alphaOlubar', 'alphaOqd1', 'alphaOqd1bar', 'alphaOqe', 'alphaOqebar', 'alphaOqq1', 'alphaOqq1bar', 'alphaOqu1', 'alphaOqu1bar', 'alphaOud1', 'alphaOud1bar', 'alphaOuH', 'alphaOuHbar', 'alphaOuu']

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
        return (-1/4 * (self.gphiVB)**(2) * (self.MVB)**(-2) + -1/4 * (self.gphiVBbar)**(2) * (self.MVB)**(-2))

    def alphaOHBt(self, ):
        return 0

    def alphaOHD(self, ):
        return (-1/2 * (self.gphiVB)**(2) * (self.MVB)**(-2) + -1 * self.gphiVB * self.gphiVBbar * (self.MVB)**(-2) + -1/2 * (self.gphiVBbar)**(2) * (self.MVB)**(-2))

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
        return -1/2 * (self.MVB)**(-2) * self.gdVB[mif1,mif2] * self.gdVB[mif3,mif4]

    def alphaOeebar(self, mif1,mif2,mif3,mif4):
        return -1/2 * (self.MVB)**(-2) * self.geVB[mif1,mif2] * self.geVB[mif3,mif4]

    def alphaOHdbar(self, mif1,mif2):
        return (-1/2 * self.gphiVB * (self.MVB)**(-2) * self.gdVB[mif1,mif2] + -1/2 * self.gphiVBbar * (self.MVB)**(-2) * self.gdVB[mif1,mif2])

    def alphaOHebar(self, mif1,mif2):
        return (-1/2 * self.gphiVB * (self.MVB)**(-2) * self.geVB[mif1,mif2] + -1/2 * self.gphiVBbar * (self.MVB)**(-2) * self.geVB[mif1,mif2])

    def alphaOHl1bar(self, mif1,mif2):
        return (-1/2 * self.gphiVB * (self.MVB)**(-2) * self.glVB[mif1,mif2] + -1/2 * self.gphiVBbar * (self.MVB)**(-2) * self.glVB[mif1,mif2])

    def alphaOHl3bar(self, mif1,mif2):
        return 0

    def alphaOHq1bar(self, mif1,mif2):
        return (-1/2 * self.gphiVB * (self.MVB)**(-2) * self.gqVB[mif1,mif2] + -1/2 * self.gphiVBbar * (self.MVB)**(-2) * self.gqVB[mif1,mif2])

    def alphaOHq3bar(self, mif1,mif2):
        return 0

    def alphaOHubar(self, mif1,mif2):
        return (-1/2 * self.gphiVB * (self.MVB)**(-2) * self.guVB[mif1,mif2] + -1/2 * self.gphiVBbar * (self.MVB)**(-2) * self.guVB[mif1,mif2])

    def alphaOllbar(self, mif1,mif2,mif3,mif4):
        return -1/2 * (self.MVB)**(-2) * self.glVB[mif1,mif2] * self.glVB[mif3,mif4]

    def alphaOuubar(self, mif1,mif2,mif3,mif4):
        return -1/2 * (self.MVB)**(-2) * self.guVB[mif1,mif2] * self.guVB[mif3,mif4]

    def alphaWeinberg(self, mif1,mif2):
        return 0

    def alphaWeinbergbar(self, mif1,mif2):
        return 0

    def alphaOdB(self, mif1,mif2):
        return 0

    def alphaOdBbar(self, mif1,mif2):
        return 0

    def alphaOdd(self, mif1,mif2,mif3,mif4):
        return -1/2 * (self.MVB)**(-2) * self.gdVB[mif1,mif2] * self.gdVB[mif3,mif4]

    def alphaOdG(self, mif1,mif2):
        return 0

    def alphaOdGbar(self, mif1,mif2):
        return 0

    def alphaOdH(self, mif1,mif2):
        return (1/4 * self.gphiVB * (self.MVB)**(-2) * self.gqVB[0,mif1] * self.yd[0,mif2] + -1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.gqVB[0,mif1] * self.yd[0,mif2] + 1/4 * self.gphiVB * (self.MVB)**(-2) * self.gqVB[mif1,0] * self.yd[0,mif2] + -1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.gqVB[mif1,0] * self.yd[0,mif2] + 1/4 * self.gphiVB * (self.MVB)**(-2) * self.gqVB[1,mif1] * self.yd[1,mif2] + -1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.gqVB[1,mif1] * self.yd[1,mif2] + 1/4 * self.gphiVB * (self.MVB)**(-2) * self.gqVB[mif1,1] * self.yd[1,mif2] + -1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.gqVB[mif1,1] * self.yd[1,mif2] + 1/4 * self.gphiVB * (self.MVB)**(-2) * self.gqVB[2,mif1] * self.yd[2,mif2] + -1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.gqVB[2,mif1] * self.yd[2,mif2] + 1/4 * self.gphiVB * (self.MVB)**(-2) * self.gqVB[mif1,2] * self.yd[2,mif2] + -1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.gqVB[mif1,2] * self.yd[2,mif2] + -1/4 * self.gphiVB * (self.MVB)**(-2) * self.gdVB[0,mif2] * self.yd[mif1,0] + 1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.gdVB[0,mif2] * self.yd[mif1,0] + -1/4 * self.gphiVB * (self.MVB)**(-2) * self.gdVB[mif2,0] * self.yd[mif1,0] + 1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.gdVB[mif2,0] * self.yd[mif1,0] + -1/4 * self.gphiVB * (self.MVB)**(-2) * self.gdVB[1,mif2] * self.yd[mif1,1] + 1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.gdVB[1,mif2] * self.yd[mif1,1] + -1/4 * self.gphiVB * (self.MVB)**(-2) * self.gdVB[mif2,1] * self.yd[mif1,1] + 1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.gdVB[mif2,1] * self.yd[mif1,1] + -1/4 * self.gphiVB * (self.MVB)**(-2) * self.gdVB[2,mif2] * self.yd[mif1,2] + 1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.gdVB[2,mif2] * self.yd[mif1,2] + -1/4 * self.gphiVB * (self.MVB)**(-2) * self.gdVB[mif2,2] * self.yd[mif1,2] + 1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.gdVB[mif2,2] * self.yd[mif1,2] + -1/4 * (self.gphiVB)**(2) * (self.MVB)**(-2) * self.yd[mif1,mif2] + 1/4 * (self.gphiVBbar)**(2) * (self.MVB)**(-2) * self.yd[mif1,mif2])

    def alphaOdHbar(self, mif1,mif2):
        return (-1/4 * self.gphiVB * (self.MVB)**(-2) * self.gqVB[0,mif1] * self.ydbar[0,mif2] + 1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.gqVB[0,mif1] * self.ydbar[0,mif2] + -1/4 * self.gphiVB * (self.MVB)**(-2) * self.gqVB[mif1,0] * self.ydbar[0,mif2] + 1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.gqVB[mif1,0] * self.ydbar[0,mif2] + -1/4 * self.gphiVB * (self.MVB)**(-2) * self.gqVB[1,mif1] * self.ydbar[1,mif2] + 1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.gqVB[1,mif1] * self.ydbar[1,mif2] + -1/4 * self.gphiVB * (self.MVB)**(-2) * self.gqVB[mif1,1] * self.ydbar[1,mif2] + 1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.gqVB[mif1,1] * self.ydbar[1,mif2] + -1/4 * self.gphiVB * (self.MVB)**(-2) * self.gqVB[2,mif1] * self.ydbar[2,mif2] + 1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.gqVB[2,mif1] * self.ydbar[2,mif2] + -1/4 * self.gphiVB * (self.MVB)**(-2) * self.gqVB[mif1,2] * self.ydbar[2,mif2] + 1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.gqVB[mif1,2] * self.ydbar[2,mif2] + 1/4 * self.gphiVB * (self.MVB)**(-2) * self.gdVB[0,mif2] * self.ydbar[mif1,0] + -1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.gdVB[0,mif2] * self.ydbar[mif1,0] + 1/4 * self.gphiVB * (self.MVB)**(-2) * self.gdVB[mif2,0] * self.ydbar[mif1,0] + -1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.gdVB[mif2,0] * self.ydbar[mif1,0] + 1/4 * self.gphiVB * (self.MVB)**(-2) * self.gdVB[1,mif2] * self.ydbar[mif1,1] + -1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.gdVB[1,mif2] * self.ydbar[mif1,1] + 1/4 * self.gphiVB * (self.MVB)**(-2) * self.gdVB[mif2,1] * self.ydbar[mif1,1] + -1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.gdVB[mif2,1] * self.ydbar[mif1,1] + 1/4 * self.gphiVB * (self.MVB)**(-2) * self.gdVB[2,mif2] * self.ydbar[mif1,2] + -1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.gdVB[2,mif2] * self.ydbar[mif1,2] + 1/4 * self.gphiVB * (self.MVB)**(-2) * self.gdVB[mif2,2] * self.ydbar[mif1,2] + -1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.gdVB[mif2,2] * self.ydbar[mif1,2] + 1/4 * (self.gphiVB)**(2) * (self.MVB)**(-2) * self.ydbar[mif1,mif2] + -1/4 * (self.gphiVBbar)**(2) * (self.MVB)**(-2) * self.ydbar[mif1,mif2])

    def alphaOdW(self, mif1,mif2):
        return 0

    def alphaOdWbar(self, mif1,mif2):
        return 0

    def alphaOeB(self, mif1,mif2):
        return 0

    def alphaOeBbar(self, mif1,mif2):
        return 0

    def alphaOed(self, mif1,mif2,mif3,mif4):
        return -1 * (self.MVB)**(-2) * self.gdVB[mif3,mif4] * self.geVB[mif1,mif2]

    def alphaOedbar(self, mif1,mif2,mif3,mif4):
        return -1 * (self.MVB)**(-2) * self.gdVB[mif3,mif4] * self.geVB[mif1,mif2]

    def alphaOee(self, mif1,mif2,mif3,mif4):
        return -1/2 * (self.MVB)**(-2) * self.geVB[mif1,mif2] * self.geVB[mif3,mif4]

    def alphaOeH(self, mif1,mif2):
        return (1/4 * self.gphiVB * (self.MVB)**(-2) * self.glVB[0,mif1] * self.yl[0,mif2] + -1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.glVB[0,mif1] * self.yl[0,mif2] + 1/4 * self.gphiVB * (self.MVB)**(-2) * self.glVB[mif1,0] * self.yl[0,mif2] + -1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.glVB[mif1,0] * self.yl[0,mif2] + 1/4 * self.gphiVB * (self.MVB)**(-2) * self.glVB[1,mif1] * self.yl[1,mif2] + -1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.glVB[1,mif1] * self.yl[1,mif2] + 1/4 * self.gphiVB * (self.MVB)**(-2) * self.glVB[mif1,1] * self.yl[1,mif2] + -1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.glVB[mif1,1] * self.yl[1,mif2] + 1/4 * self.gphiVB * (self.MVB)**(-2) * self.glVB[2,mif1] * self.yl[2,mif2] + -1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.glVB[2,mif1] * self.yl[2,mif2] + 1/4 * self.gphiVB * (self.MVB)**(-2) * self.glVB[mif1,2] * self.yl[2,mif2] + -1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.glVB[mif1,2] * self.yl[2,mif2] + -1/4 * self.gphiVB * (self.MVB)**(-2) * self.geVB[0,mif2] * self.yl[mif1,0] + 1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.geVB[0,mif2] * self.yl[mif1,0] + -1/4 * self.gphiVB * (self.MVB)**(-2) * self.geVB[mif2,0] * self.yl[mif1,0] + 1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.geVB[mif2,0] * self.yl[mif1,0] + -1/4 * self.gphiVB * (self.MVB)**(-2) * self.geVB[1,mif2] * self.yl[mif1,1] + 1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.geVB[1,mif2] * self.yl[mif1,1] + -1/4 * self.gphiVB * (self.MVB)**(-2) * self.geVB[mif2,1] * self.yl[mif1,1] + 1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.geVB[mif2,1] * self.yl[mif1,1] + -1/4 * self.gphiVB * (self.MVB)**(-2) * self.geVB[2,mif2] * self.yl[mif1,2] + 1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.geVB[2,mif2] * self.yl[mif1,2] + -1/4 * self.gphiVB * (self.MVB)**(-2) * self.geVB[mif2,2] * self.yl[mif1,2] + 1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.geVB[mif2,2] * self.yl[mif1,2] + -1/4 * (self.gphiVB)**(2) * (self.MVB)**(-2) * self.yl[mif1,mif2] + 1/4 * (self.gphiVBbar)**(2) * (self.MVB)**(-2) * self.yl[mif1,mif2])

    def alphaOeHbar(self, mif1,mif2):
        return (-1/4 * self.gphiVB * (self.MVB)**(-2) * self.glVB[0,mif1] * self.ylbar[0,mif2] + 1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.glVB[0,mif1] * self.ylbar[0,mif2] + -1/4 * self.gphiVB * (self.MVB)**(-2) * self.glVB[mif1,0] * self.ylbar[0,mif2] + 1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.glVB[mif1,0] * self.ylbar[0,mif2] + -1/4 * self.gphiVB * (self.MVB)**(-2) * self.glVB[1,mif1] * self.ylbar[1,mif2] + 1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.glVB[1,mif1] * self.ylbar[1,mif2] + -1/4 * self.gphiVB * (self.MVB)**(-2) * self.glVB[mif1,1] * self.ylbar[1,mif2] + 1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.glVB[mif1,1] * self.ylbar[1,mif2] + -1/4 * self.gphiVB * (self.MVB)**(-2) * self.glVB[2,mif1] * self.ylbar[2,mif2] + 1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.glVB[2,mif1] * self.ylbar[2,mif2] + -1/4 * self.gphiVB * (self.MVB)**(-2) * self.glVB[mif1,2] * self.ylbar[2,mif2] + 1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.glVB[mif1,2] * self.ylbar[2,mif2] + 1/4 * self.gphiVB * (self.MVB)**(-2) * self.geVB[0,mif2] * self.ylbar[mif1,0] + -1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.geVB[0,mif2] * self.ylbar[mif1,0] + 1/4 * self.gphiVB * (self.MVB)**(-2) * self.geVB[mif2,0] * self.ylbar[mif1,0] + -1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.geVB[mif2,0] * self.ylbar[mif1,0] + 1/4 * self.gphiVB * (self.MVB)**(-2) * self.geVB[1,mif2] * self.ylbar[mif1,1] + -1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.geVB[1,mif2] * self.ylbar[mif1,1] + 1/4 * self.gphiVB * (self.MVB)**(-2) * self.geVB[mif2,1] * self.ylbar[mif1,1] + -1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.geVB[mif2,1] * self.ylbar[mif1,1] + 1/4 * self.gphiVB * (self.MVB)**(-2) * self.geVB[2,mif2] * self.ylbar[mif1,2] + -1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.geVB[2,mif2] * self.ylbar[mif1,2] + 1/4 * self.gphiVB * (self.MVB)**(-2) * self.geVB[mif2,2] * self.ylbar[mif1,2] + -1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.geVB[mif2,2] * self.ylbar[mif1,2] + 1/4 * (self.gphiVB)**(2) * (self.MVB)**(-2) * self.ylbar[mif1,mif2] + -1/4 * (self.gphiVBbar)**(2) * (self.MVB)**(-2) * self.ylbar[mif1,mif2])

    def alphaOeu(self, mif1,mif2,mif3,mif4):
        return -1 * (self.MVB)**(-2) * self.geVB[mif1,mif2] * self.guVB[mif3,mif4]

    def alphaOeubar(self, mif1,mif2,mif3,mif4):
        return -1 * (self.MVB)**(-2) * self.geVB[mif1,mif2] * self.guVB[mif3,mif4]

    def alphaOeW(self, mif1,mif2):
        return 0

    def alphaOeWbar(self, mif1,mif2):
        return 0

    def alphaOHd(self, mif1,mif2):
        return (-1/2 * self.gphiVB * (self.MVB)**(-2) * self.gdVB[mif1,mif2] + -1/2 * self.gphiVBbar * (self.MVB)**(-2) * self.gdVB[mif1,mif2])

    def alphaOHe(self, mif1,mif2):
        return (-1/2 * self.gphiVB * (self.MVB)**(-2) * self.geVB[mif1,mif2] + -1/2 * self.gphiVBbar * (self.MVB)**(-2) * self.geVB[mif1,mif2])

    def alphaOHl1(self, mif1,mif2):
        return (-1/2 * self.gphiVB * (self.MVB)**(-2) * self.glVB[mif1,mif2] + -1/2 * self.gphiVBbar * (self.MVB)**(-2) * self.glVB[mif1,mif2])

    def alphaOHl3(self, mif1,mif2):
        return 0

    def alphaOHq1(self, mif1,mif2):
        return (-1/2 * self.gphiVB * (self.MVB)**(-2) * self.gqVB[mif1,mif2] + -1/2 * self.gphiVBbar * (self.MVB)**(-2) * self.gqVB[mif1,mif2])

    def alphaOHq3(self, mif1,mif2):
        return 0

    def alphaOHu(self, mif1,mif2):
        return (-1/2 * self.gphiVB * (self.MVB)**(-2) * self.guVB[mif1,mif2] + -1/2 * self.gphiVBbar * (self.MVB)**(-2) * self.guVB[mif1,mif2])

    def alphaOHud(self, mif1,mif2):
        return 0

    def alphaOHudbar(self, mif1,mif2):
        return 0

    def alphaOld(self, mif1,mif2,mif3,mif4):
        return -1 * (self.MVB)**(-2) * self.gdVB[mif3,mif4] * self.glVB[mif1,mif2]

    def alphaOldbar(self, mif1,mif2,mif3,mif4):
        return -1 * (self.MVB)**(-2) * self.gdVB[mif3,mif4] * self.glVB[mif1,mif2]

    def alphaOle(self, mif1,mif2,mif3,mif4):
        return -1 * (self.MVB)**(-2) * self.geVB[mif3,mif4] * self.glVB[mif1,mif2]

    def alphaOlebar(self, mif1,mif2,mif3,mif4):
        return -1 * (self.MVB)**(-2) * self.geVB[mif3,mif4] * self.glVB[mif1,mif2]

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
        return -1/2 * (self.MVB)**(-2) * self.glVB[mif1,mif2] * self.glVB[mif3,mif4]

    def alphaOlq1(self, mif1,mif2,mif3,mif4):
        return -1 * (self.MVB)**(-2) * self.glVB[mif1,mif2] * self.gqVB[mif3,mif4]

    def alphaOlq1bar(self, mif1,mif2,mif3,mif4):
        return -1 * (self.MVB)**(-2) * self.glVB[mif1,mif2] * self.gqVB[mif3,mif4]

    def alphaOlq3(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOlq3bar(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOlu(self, mif1,mif2,mif3,mif4):
        return -1 * (self.MVB)**(-2) * self.glVB[mif1,mif2] * self.guVB[mif3,mif4]

    def alphaOlubar(self, mif1,mif2,mif3,mif4):
        return -1 * (self.MVB)**(-2) * self.glVB[mif1,mif2] * self.guVB[mif3,mif4]

    def alphaOqd1(self, mif1,mif2,mif3,mif4):
        return -1 * (self.MVB)**(-2) * self.gdVB[mif3,mif4] * self.gqVB[mif1,mif2]

    def alphaOqd1bar(self, mif1,mif2,mif3,mif4):
        return -1 * (self.MVB)**(-2) * self.gdVB[mif3,mif4] * self.gqVB[mif1,mif2]

    def alphaOqd8(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOqd8bar(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOqe(self, mif1,mif2,mif3,mif4):
        return -1 * (self.MVB)**(-2) * self.geVB[mif3,mif4] * self.gqVB[mif1,mif2]

    def alphaOqebar(self, mif1,mif2,mif3,mif4):
        return -1 * (self.MVB)**(-2) * self.geVB[mif3,mif4] * self.gqVB[mif1,mif2]

    def alphaOqq1(self, mif1,mif2,mif3,mif4):
        return -1/2 * (self.MVB)**(-2) * self.gqVB[mif1,mif2] * self.gqVB[mif3,mif4]

    def alphaOqq1bar(self, mif1,mif2,mif3,mif4):
        return -1/2 * (self.MVB)**(-2) * self.gqVB[mif1,mif2] * self.gqVB[mif3,mif4]

    def alphaOqq3(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOqq3bar(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOqu1(self, mif1,mif2,mif3,mif4):
        return -1 * (self.MVB)**(-2) * self.gqVB[mif1,mif2] * self.guVB[mif3,mif4]

    def alphaOqu1bar(self, mif1,mif2,mif3,mif4):
        return -1 * (self.MVB)**(-2) * self.gqVB[mif1,mif2] * self.guVB[mif3,mif4]

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
        return -1 * (self.MVB)**(-2) * self.gdVB[mif3,mif4] * self.guVB[mif1,mif2]

    def alphaOud1bar(self, mif1,mif2,mif3,mif4):
        return -1 * (self.MVB)**(-2) * self.gdVB[mif3,mif4] * self.guVB[mif1,mif2]

    def alphaOud8(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOud8bar(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOuG(self, mif1,mif2):
        return 0

    def alphaOuGbar(self, mif1,mif2):
        return 0

    def alphaOuH(self, mif1,mif2):
        return (1/4 * self.gphiVB * (self.MVB)**(-2) * self.gqVB[0,mif1] * self.yu[0,mif2] + -1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.gqVB[0,mif1] * self.yu[0,mif2] + 1/4 * self.gphiVB * (self.MVB)**(-2) * self.gqVB[mif1,0] * self.yu[0,mif2] + -1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.gqVB[mif1,0] * self.yu[0,mif2] + 1/4 * self.gphiVB * (self.MVB)**(-2) * self.gqVB[1,mif1] * self.yu[1,mif2] + -1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.gqVB[1,mif1] * self.yu[1,mif2] + 1/4 * self.gphiVB * (self.MVB)**(-2) * self.gqVB[mif1,1] * self.yu[1,mif2] + -1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.gqVB[mif1,1] * self.yu[1,mif2] + 1/4 * self.gphiVB * (self.MVB)**(-2) * self.gqVB[2,mif1] * self.yu[2,mif2] + -1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.gqVB[2,mif1] * self.yu[2,mif2] + 1/4 * self.gphiVB * (self.MVB)**(-2) * self.gqVB[mif1,2] * self.yu[2,mif2] + -1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.gqVB[mif1,2] * self.yu[2,mif2] + -1/4 * self.gphiVB * (self.MVB)**(-2) * self.guVB[0,mif2] * self.yu[mif1,0] + 1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.guVB[0,mif2] * self.yu[mif1,0] + -1/4 * self.gphiVB * (self.MVB)**(-2) * self.guVB[mif2,0] * self.yu[mif1,0] + 1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.guVB[mif2,0] * self.yu[mif1,0] + -1/4 * self.gphiVB * (self.MVB)**(-2) * self.guVB[1,mif2] * self.yu[mif1,1] + 1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.guVB[1,mif2] * self.yu[mif1,1] + -1/4 * self.gphiVB * (self.MVB)**(-2) * self.guVB[mif2,1] * self.yu[mif1,1] + 1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.guVB[mif2,1] * self.yu[mif1,1] + -1/4 * self.gphiVB * (self.MVB)**(-2) * self.guVB[2,mif2] * self.yu[mif1,2] + 1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.guVB[2,mif2] * self.yu[mif1,2] + -1/4 * self.gphiVB * (self.MVB)**(-2) * self.guVB[mif2,2] * self.yu[mif1,2] + 1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.guVB[mif2,2] * self.yu[mif1,2] + 1/4 * (self.gphiVB)**(2) * (self.MVB)**(-2) * self.yu[mif1,mif2] + -1/4 * (self.gphiVBbar)**(2) * (self.MVB)**(-2) * self.yu[mif1,mif2])

    def alphaOuHbar(self, mif1,mif2):
        return (-1/4 * self.gphiVB * (self.MVB)**(-2) * self.gqVB[0,mif1] * self.yubar[0,mif2] + 1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.gqVB[0,mif1] * self.yubar[0,mif2] + -1/4 * self.gphiVB * (self.MVB)**(-2) * self.gqVB[mif1,0] * self.yubar[0,mif2] + 1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.gqVB[mif1,0] * self.yubar[0,mif2] + -1/4 * self.gphiVB * (self.MVB)**(-2) * self.gqVB[1,mif1] * self.yubar[1,mif2] + 1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.gqVB[1,mif1] * self.yubar[1,mif2] + -1/4 * self.gphiVB * (self.MVB)**(-2) * self.gqVB[mif1,1] * self.yubar[1,mif2] + 1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.gqVB[mif1,1] * self.yubar[1,mif2] + -1/4 * self.gphiVB * (self.MVB)**(-2) * self.gqVB[2,mif1] * self.yubar[2,mif2] + 1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.gqVB[2,mif1] * self.yubar[2,mif2] + -1/4 * self.gphiVB * (self.MVB)**(-2) * self.gqVB[mif1,2] * self.yubar[2,mif2] + 1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.gqVB[mif1,2] * self.yubar[2,mif2] + 1/4 * self.gphiVB * (self.MVB)**(-2) * self.guVB[0,mif2] * self.yubar[mif1,0] + -1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.guVB[0,mif2] * self.yubar[mif1,0] + 1/4 * self.gphiVB * (self.MVB)**(-2) * self.guVB[mif2,0] * self.yubar[mif1,0] + -1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.guVB[mif2,0] * self.yubar[mif1,0] + 1/4 * self.gphiVB * (self.MVB)**(-2) * self.guVB[1,mif2] * self.yubar[mif1,1] + -1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.guVB[1,mif2] * self.yubar[mif1,1] + 1/4 * self.gphiVB * (self.MVB)**(-2) * self.guVB[mif2,1] * self.yubar[mif1,1] + -1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.guVB[mif2,1] * self.yubar[mif1,1] + 1/4 * self.gphiVB * (self.MVB)**(-2) * self.guVB[2,mif2] * self.yubar[mif1,2] + -1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.guVB[2,mif2] * self.yubar[mif1,2] + 1/4 * self.gphiVB * (self.MVB)**(-2) * self.guVB[mif2,2] * self.yubar[mif1,2] + -1/4 * self.gphiVBbar * (self.MVB)**(-2) * self.guVB[mif2,2] * self.yubar[mif1,2] + -1/4 * (self.gphiVB)**(2) * (self.MVB)**(-2) * self.yubar[mif1,mif2] + 1/4 * (self.gphiVBbar)**(2) * (self.MVB)**(-2) * self.yubar[mif1,mif2])

    def alphaOuu(self, mif1,mif2,mif3,mif4):
        return -1/2 * (self.MVB)**(-2) * self.guVB[mif1,mif2] * self.guVB[mif3,mif4]

    def alphaOuW(self, mif1,mif2):
        return 0

    def alphaOuWbar(self, mif1,mif2):
        return 0
