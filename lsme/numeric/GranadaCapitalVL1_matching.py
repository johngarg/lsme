## Symbols used:
# ydbar
# yd
# yl
# ylbar
# muH
# lam
# yubar
# yu

## 2025-10-08 JG: Copied by hand from 1711.10391v2

import sys
import os
import numpy as np

import lsme.numeric.matchingresult as matchingresult

class GranadaVL1MatchingResult(matchingresult.GenericMatchingResult):
    def __init__(self, name='VL1'):
        super().__init__(name)
        self.MVL1 = 1
        self.gammaVL1 = 1
        self.gammaVL1bar = 1
        self.gBVL1 = 1
        self.gBVL1bar = 1
        self.gWVL1 = 1
        self.gWVL1bar = 1
        self.gBtildeVL1 = 1
        self.gBtildeVL1bar = 1
        self.gWtildeVL1 = 1
        self.gWtildeVL1bar = 1
        self.h1VL1 = 1
        self.h1VL1bar = 1
        self.h2VL1 = 1
        self.h2VL1bar = 1
        self.h3VL1 = 1
        self.h3VL1bar = 1
        self.exotic_params = ['MVL1', 'gammaVL1', 'gammaVL1bar', 'gBVL1', 'gBVL1bar', 'gWVL1', 'gWVL1bar', 'gBtildeVL1', 'gBtildeVL1bar', 'gWtildeVL1', 'gWtildeVL1bar', 'h1VL1', 'h1VL1bar', 'h2VL1', 'h2VL1bar', 'h3VL1', 'h3VL1bar']
        self.nonvanishing = ['alphaOH', 'alphaOHB', 'alphaOHBox', 'alphaOHBt', 'alphaOHD', 'alphaOHW', 'alphaOHWB', 'alphaOHWBt', 'alphaOHWt', 'alphaOdH', 'alphaOdHbar', 'alphaOeH', 'alphaOeHbar' 'alphaOHd', 'alphaOHe', 'alphaOHl1', 'alphaOHl3', 'alphaOHq1', 'alphaOHq3', 'alphaOuH', 'alphaOuHbar']

    def alphaO3G(self, ):
        return 0

    def alphaO3Gt(self, ):
        return 0

    def alphaO3W(self, ):
        return 0

    def alphaO3Wt(self, ):
        return 0

    def alphaOH(self, ):
        term_1 = 2 * self.g2 * self.lam * self.gWVL1 * self.gammaVL1 * self.gammaVL1bar * (self.MVL1)**(-4)
        term_2 = - 2 * self.lam * self.h1VL1 * self.gammaVL1 * self.gammaVL1bar * (self.MVL1)**(-4)
        return term_1 + term_2

    def alphaOHB(self, ):
        term_1 = - self.g1**2 * self.gammaVL1 * self.gammaVL1bar * (self.MVL1)**(-4) / 8.0
        term_2 = - self.g1 * self.gBVL1 *  self.gammaVL1 * self.gammaVL1bar * (self.MVL1)**(-4) / 4.0
        return term_1 + term_2

    def alphaOHBox(self, ):
        term_1 = self.g1 * self.gBVL1 * self.gammaVL1 * self.gammaVL1bar * (self.MVL1)**(-4) / 4.0
        term_2 = 3 * self.g2 * self.gWVL1 * self.gammaVL1 * self.gammaVL1bar * (self.MVL1)**(-4) / 4.0
        term_3 = - self.h1VL1 * self.gammaVL1 * self.gammaVL1bar * (self.MVL1)**(-4) / 2.0
        term_4 = np.real(self.h3VL1 * self.gammaVL1bar * self.gammaVL1bar) * (self.MVL1)**(-4)
        return term_1 + term_2 + term_3 + term_4

    def alphaOHBt(self, ):
        return - self.g1 * self.gBtildeVL1 *  self.gammaVL1 * self.gammaVL1bar * (self.MVL1)**(-4) / 4.0

    def alphaOHD(self, ):
        term_1 = self.g1 * self.gBVL1 * self.gammaVL1 * self.gammaVL1bar * (self.MVL1)**(-4)
        term_2 = - self.h2VL1 * self.gammaVL1 * self.gammaVL1bar * (self.MVL1)**(-4)
        term_3 = 2 * np.real(self.h3VL1 * self.gammaVL1bar * self.gammaVL1bar) * (self.MVL1)**(-4)
        return term_1 + term_2 + term_3

    def alphaOHG(self, ):
        return 0

    def alphaOHGt(self, ):
        return 0

    def alphaOHW(self, ):
        term_1 = - self.g2**2 * self.gammaVL1 * self.gammaVL1bar * (self.MVL1)**(-4) / 8.0
        term_2 = - self.g2 * self.gWVL1 *  self.gammaVL1 * self.gammaVL1bar * (self.MVL1)**(-4) / 4.0
        return term_1 + term_2

    def alphaOHWB(self, ):
        term_1 = - self.g1 * self.g2 * self.gammaVL1 * self.gammaVL1bar * (self.MVL1)**(-4) / 4.0
        term_2 = - self.g2 * self.gBVL1 * self.gammaVL1 * self.gammaVL1bar * (self.MVL1)**(-4) / 4.0
        term_3 = - self.g1 * self.gWVL1 * self.gammaVL1 * self.gammaVL1bar * (self.MVL1)**(-4) / 4.0
        return term_1 + term_2 + term_3

    def alphaOHWBt(self, ):
        term_1 = - self.g2 * self.gBtildeVL1 * self.gammaVL1 * self.gammaVL1bar * (self.MVL1)**(-4) / 4.0
        term_2 = - self.g1 * self.gWtildeVL1 * self.gammaVL1 * self.gammaVL1bar * (self.MVL1)**(-4) / 4.0
        return term_1 + term_2

    def alphaOHWt(self, ):
        return - self.g2 * self.gWtildeVL1 *  self.gammaVL1 * self.gammaVL1bar * (self.MVL1)**(-4) / 4.0

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
        term_1 = self.g2 * self.gWVL1 * self.gammaVL1 * self.gammaVL1bar * (self.MVL1)**(-4) / 2.0
        term_2 = - self.h1VL1 * self.gammaVL1 * self.gammaVL1bar * (self.MVL1)**(-4) / 2.0
        term_3 = - 1j * np.imag(self.h3VL1 * self.gammaVL1bar * self.gammaVL1bar) * (self.MVL1)**(-4)
        return np.conj(self.yd[mif2,mif1]) * (term_1 + term_2 + term_3)

    def alphaOdHbar(self, mif1,mif2):
        return np.conj(self.alphaOdH(mif1, mif2))

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
        term_1 = self.g2 * self.gWVL1 * self.gammaVL1 * self.gammaVL1bar * (self.MVL1)**(-4) / 2.0
        term_2 = - self.h1VL1 * self.gammaVL1 * self.gammaVL1bar * (self.MVL1)**(-4) / 2.0
        term_3 = - 1j * np.imag(self.h3VL1 * self.gammaVL1bar * self.gammaVL1bar) * (self.MVL1)**(-4)
        return np.conj(self.yl[mif2,mif1]) * (term_1 + term_2 + term_3)

    def alphaOeHbar(self, mif1,mif2):
        return np.conj(self.alphaOeH(mif1, mif2))

    def alphaOeu(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOeubar(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOeW(self, mif1,mif2):
        return 0

    def alphaOeWbar(self, mif1,mif2):
        return 0

    def alphaOHd(self, mif1,mif2):
        return - self.g1 * np.kron(mif1,mif2) * self.gBVL1 * self.gammaVL1 * self.gammaVL1bar * (self.MVL1)**(-4) / 6.0

    def alphaOHe(self, mif1,mif2):
        return - self.g1 * np.kron(mif1,mif2) * self.gBVL1 * self.gammaVL1 * self.gammaVL1bar * (self.MVL1)**(-4) / 2.0

    def alphaOHl1(self, mif1,mif2):
        return - self.g1 * np.kron(mif1,mif2) * self.gBVL1 * self.gammaVL1 * self.gammaVL1bar * (self.MVL1)**(-4) / 4.0

    def alphaOHl3(self, mif1,mif2):
        return self.g2 * np.kron(mif1,mif2) * self.gWVL1 * self.gammaVL1 * self.gammaVL1bar * (self.MVL1)**(-4) / 4.0

    def alphaOHq1(self, mif1,mif2):
        return self.g1 * np.kron(mif1,mif2) * self.gBVL1 * self.gammaVL1 * self.gammaVL1bar * (self.MVL1)**(-4) / 12.0

    def alphaOHq3(self, mif1,mif2):
        return self.g2 * np.kron(mif1,mif2) * self.gWVL1 * self.gammaVL1 * self.gammaVL1bar * (self.MVL1)**(-4) / 4.0

    def alphaOHu(self, mif1,mif2):
        return self.g1 * np.kron(mif1,mif2) * self.gBVL1 * self.gammaVL1 * self.gammaVL1bar * (self.MVL1)**(-4) / 3.0

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
        return 0

    def alphaOqq3bar(self, mif1,mif2,mif3,mif4):
        return 0

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
        term_1 = self.g2 * self.gWVL1 * self.gammaVL1 * self.gammaVL1bar * (self.MVL1)**(-4) / 2.0
        term_2 = - self.h1VL1 * self.gammaVL1 * self.gammaVL1bar * (self.MVL1)**(-4) / 2.0
        term_3 = - 1j * np.imag(self.h3VL1 * self.gammaVL1bar * self.gammaVL1bar) * (self.MVL1)**(-4)
        return np.conj(self.yu[mif2,mif1]) * np.conj(term_1 + term_2 + term_3)

    def alphaOuHbar(self, mif1,mif2):
        return np.conj(self.alphaOuH(mif1,mif2))

    def alphaOuu(self, mif1,mif2,mif3,mif4):
        return 0

    def alphaOuW(self, mif1,mif2):
        return 0

    def alphaOuWbar(self, mif1,mif2):
        return 0
