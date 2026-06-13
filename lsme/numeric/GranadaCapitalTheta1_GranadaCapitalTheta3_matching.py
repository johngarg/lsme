## Generated from SOLD output Attachments/theta1theta3.wl.
## Coupling normalisations follow 2412.01759 Table 2 and 20260612T155918.

import sys
import os
import numpy as np

import lsme.numeric.matchingresult as matchingresult


class GranadaTheta1_GranadaTheta3MatchingResult(matchingresult.GenericMatchingResult):
    def __init__(self, name='Theta1_Theta3'):
        super().__init__(name)
        self.MTheta1 = 1
        self.MTheta3 = 1
        self.lambdaTheta1 = 1
        self.lambdaTheta1bar = 1
        self.lambdaTheta3 = 1
        self.lambdaTheta3bar = 1
        self.lambdaHatTheta1 = 1
        self.lambdaHatTheta1bar = 1
        self.lambdaHatTheta3 = 1
        self.lambdaHatTheta3bar = 1
        self.lambdaHatPrimeTheta1 = 1
        self.lambdaHatPrimeTheta1bar = 1
        self.lambdaHatPrimeTheta3 = 1
        self.lambdaHatPrimeTheta3bar = 1
        self.lambdaHatPrimePrimeTheta1 = 1
        self.lambdaHatPrimePrimeTheta1bar = 1
        self.lambdaHatTheta1Theta3 = 1
        self.lambdaHatTheta1Theta3bar = 1
        self.exotic_params = ['MTheta1', 'MTheta3', 'lambdaTheta1', 'lambdaTheta1bar', 'lambdaTheta3', 'lambdaTheta3bar', 'lambdaHatTheta1', 'lambdaHatTheta1bar', 'lambdaHatTheta3', 'lambdaHatTheta3bar', 'lambdaHatPrimeTheta1', 'lambdaHatPrimeTheta1bar', 'lambdaHatPrimeTheta3', 'lambdaHatPrimeTheta3bar', 'lambdaHatPrimePrimeTheta1', 'lambdaHatPrimePrimeTheta1bar', 'lambdaHatTheta1Theta3', 'lambdaHatTheta1Theta3bar']
        self.nonvanishing = ['alphaO3W', 'alphaOHW', 'alphaOHB', 'alphaOHWB', 'alphaOHBox', 'alphaOHD', 'alphaOHq1', 'alphaOHq3', 'alphaOHu', 'alphaOHd', 'alphaOHl1', 'alphaOHl3', 'alphaOHe', 'alphaOuH', 'alphaOdH', 'alphaOeH', 'alphaOqq1', 'alphaOqq3', 'alphaOuu', 'alphaOdd', 'alphaOud1', 'alphaOqu1', 'alphaOqd1', 'alphaOll', 'alphaOee', 'alphaOle', 'alphaOlq1', 'alphaOlq3', 'alphaOeu', 'alphaOed', 'alphaOqe', 'alphaOlu', 'alphaOld']

    def alphaO3G(self, ):
        return 0

    def alphaO3Gt(self, ):
        return 0

    def alphaO3W(self, ):
        return (1/576 * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) * (self.g2)**(3) + 1/576 * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(-2) * (self.g2)**(3))

    def alphaO3Wt(self, ):
        return 0

    def alphaOHG(self, ):
        return 0

    def alphaOHGt(self, ):
        return 0

    def alphaOHW(self, ):
        return (-5/192 * self.lambdaHatTheta1 * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) * (self.g2)**(2) + -5/192 * self.lambdaHatTheta3 * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(-2) * (self.g2)**(2))

    def alphaOHWt(self, ):
        return 0

    def alphaOHB(self, ):
        return (-3/64 * self.lambdaHatTheta3 * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(-2) * (self.g1)**(2) + -1/192 * self.lambdaHatTheta1 * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) * (self.g1)**(2))

    def alphaOHBt(self, ):
        return 0

    def alphaOHWB(self, ):
        return (-5/64 * self.g1 * self.g2 * self.lambdaHatPrimeTheta3 * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(-2) + -5/192 * self.g1 * self.g2 * self.lambdaHatPrimeTheta1 * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2))

    def alphaOHWBt(self, ):
        return 0

    def alphaOHBox(self, ):
        return (-3/1280 * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(-2) * (self.g1)**(4) + -1/48 * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) * (self.lambdaHatTheta1)**(2) + -1/48 * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(-2) * (self.lambdaHatTheta3)**(2) + -1/256 * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) * (self.g2)**(4) + -1/256 * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(-2) * (self.g2)**(4) + -1/3840 * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) * (self.g1)**(4) + 5/192 * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) * (self.lambdaHatPrimeTheta1)**(2) + 5/192 * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(-2) * (self.lambdaHatPrimeTheta3)**(2) + 1/16 * self.lambdaTheta1 * self.lambdaTheta1bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) + 1/18 * self.lambdaHatPrimePrimeTheta1 * self.lambdaHatPrimePrimeTheta1bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) + 3/16 * self.lambdaTheta3 * self.lambdaTheta3bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(-2) + -15/16 * self.lambdaHatTheta1Theta3 * self.lambdaHatTheta1Theta3bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(4) * (((self.MTheta1)**(2) + -1 * (self.MTheta3)**(2)))**(-3) + -5/8 * self.lambdaHatTheta1Theta3 * self.lambdaHatTheta1Theta3bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(2) * (((self.MTheta1)**(2) + -1 * (self.MTheta3)**(2)))**(-2) + -5/8 * self.lambdaHatTheta1Theta3 * self.lambdaHatTheta1Theta3bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(2) * (((self.MTheta1)**(2) + -1 * (self.MTheta3)**(2)))**(-2) + 15/16 * self.lambdaHatTheta1Theta3 * self.lambdaHatTheta1Theta3bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(4) * (((self.MTheta1)**(2) + -1 * (self.MTheta3)**(2)))**(-3) + -5/8 * self.lambdaHatTheta1Theta3 * self.lambdaHatTheta1Theta3bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(4) * (((self.MTheta1)**(2) + -1 * (self.MTheta3)**(2)))**(-3) * np.log((self.MTheta1)**(2) * (self.mu)**(-2)) + 5/8 * self.lambdaHatTheta1Theta3 * self.lambdaHatTheta1Theta3bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(2) * (((self.MTheta1)**(2) + -1 * (self.MTheta3)**(2)))**(-2) * np.log((self.MTheta1)**(2) * (self.mu)**(-2)) + 5/8 * self.lambdaHatTheta1Theta3 * self.lambdaHatTheta1Theta3bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(2) * (((self.MTheta1)**(2) + -1 * (self.MTheta3)**(2)))**(-2) * np.log((self.MTheta3)**(2) * (self.mu)**(-2)) + 5/8 * self.lambdaHatTheta1Theta3 * self.lambdaHatTheta1Theta3bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(4) * (((self.MTheta1)**(2) + -1 * (self.MTheta3)**(2)))**(-3) * np.log((self.MTheta3)**(2) * (self.mu)**(-2)))

    def alphaOHD(self, ):
        return (-5/48 * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) * (self.lambdaHatPrimeTheta1)**(2) + -5/48 * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(-2) * (self.lambdaHatPrimeTheta3)**(2) + -3/320 * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(-2) * (self.g1)**(4) + -1/960 * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) * (self.g1)**(4) + -1/8 * self.lambdaTheta1 * self.lambdaTheta1bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) + 1/9 * self.lambdaHatPrimePrimeTheta1 * self.lambdaHatPrimePrimeTheta1bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) + 3/8 * self.lambdaTheta3 * self.lambdaTheta3bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(-2) + -15/8 * self.lambdaHatTheta1Theta3 * self.lambdaHatTheta1Theta3bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(4) * (((self.MTheta1)**(2) + -1 * (self.MTheta3)**(2)))**(-3) + -5/4 * self.lambdaHatTheta1Theta3 * self.lambdaHatTheta1Theta3bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(2) * (((self.MTheta1)**(2) + -1 * (self.MTheta3)**(2)))**(-2) + -5/4 * self.lambdaHatTheta1Theta3 * self.lambdaHatTheta1Theta3bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(2) * (((self.MTheta1)**(2) + -1 * (self.MTheta3)**(2)))**(-2) + 15/8 * self.lambdaHatTheta1Theta3 * self.lambdaHatTheta1Theta3bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(4) * (((self.MTheta1)**(2) + -1 * (self.MTheta3)**(2)))**(-3) + -5/4 * self.lambdaHatTheta1Theta3 * self.lambdaHatTheta1Theta3bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(4) * (((self.MTheta1)**(2) + -1 * (self.MTheta3)**(2)))**(-3) * np.log((self.MTheta1)**(2) * (self.mu)**(-2)) + 5/4 * self.lambdaHatTheta1Theta3 * self.lambdaHatTheta1Theta3bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(2) * (((self.MTheta1)**(2) + -1 * (self.MTheta3)**(2)))**(-2) * np.log((self.MTheta1)**(2) * (self.mu)**(-2)) + 5/4 * self.lambdaHatTheta1Theta3 * self.lambdaHatTheta1Theta3bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(2) * (((self.MTheta1)**(2) + -1 * (self.MTheta3)**(2)))**(-2) * np.log((self.MTheta3)**(2) * (self.mu)**(-2)) + 5/4 * self.lambdaHatTheta1Theta3 * self.lambdaHatTheta1Theta3bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(4) * (((self.MTheta1)**(2) + -1 * (self.MTheta3)**(2)))**(-3) * np.log((self.MTheta3)**(2) * (self.mu)**(-2)))

    def alphaOH(self, ):
        return 0

    def alphaOuG(self, ):
        return 0

    def alphaOuW(self, ):
        return 0

    def alphaOuB(self, ):
        return 0

    def alphaOdG(self, ):
        return 0

    def alphaOdW(self, ):
        return 0

    def alphaOdB(self, ):
        return 0

    def alphaOeW(self, ):
        return 0

    def alphaOeB(self, ):
        return 0

    def alphaOHq1(self, fla,flb):
        return (-1/640 * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(-2) * (self.g1)**(4) * self.kronecker_delta(fla,flb) + -1/5760 * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) * (self.g1)**(4) * self.kronecker_delta(fla,flb))

    def alphaOHq3(self, fla,flb):
        return (-1/384 * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) * (self.g2)**(4) * self.kronecker_delta(fla,flb) + -1/384 * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(-2) * (self.g2)**(4) * self.kronecker_delta(fla,flb))

    def alphaOHu(self, fla,flb):
        return (-1/160 * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(-2) * (self.g1)**(4) * self.kronecker_delta(fla,flb) + -1/1440 * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) * (self.g1)**(4) * self.kronecker_delta(fla,flb))

    def alphaOHd(self, fla,flb):
        return (1/320 * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(-2) * (self.g1)**(4) * self.kronecker_delta(fla,flb) + 1/2880 * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) * (self.g1)**(4) * self.kronecker_delta(fla,flb))

    def alphaOHud(self, ):
        return 0

    def alphaOHl1(self, fla,flb):
        return (1/1920 * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) * (self.g1)**(4) * self.kronecker_delta(fla,flb) + 3/640 * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(-2) * (self.g1)**(4) * self.kronecker_delta(fla,flb))

    def alphaOHl3(self, fla,flb):
        return (-1/384 * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) * (self.g2)**(4) * self.kronecker_delta(fla,flb) + -1/384 * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(-2) * (self.g2)**(4) * self.kronecker_delta(fla,flb))

    def alphaOHe(self, fla,flb):
        return (1/960 * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) * (self.g1)**(4) * self.kronecker_delta(fla,flb) + 3/320 * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(-2) * (self.g1)**(4) * self.kronecker_delta(fla,flb))

    def alphaOuH(self, fla,flb):
        return (-1/384 * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) * (self.g2)**(4) * self.yu[fla,flb] + -1/384 * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(-2) * (self.g2)**(4) * self.yu[fla,flb] + 5/96 * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) * (self.lambdaHatPrimeTheta1)**(2) * self.yu[fla,flb] + 5/96 * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(-2) * (self.lambdaHatPrimeTheta3)**(2) * self.yu[fla,flb] + 1/18 * self.lambdaHatPrimePrimeTheta1 * self.lambdaHatPrimePrimeTheta1bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) * self.yu[fla,flb] + 3/16 * self.lambdaTheta3 * self.lambdaTheta3bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(-2) * self.yu[fla,flb] + 5/48 * self.lambdaTheta1 * self.lambdaTheta1bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) * self.yu[fla,flb] + -15/16 * self.lambdaHatTheta1Theta3 * self.lambdaHatTheta1Theta3bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(4) * (((self.MTheta1)**(2) + -1 * (self.MTheta3)**(2)))**(-3) * self.yu[fla,flb] + -5/8 * self.lambdaHatTheta1Theta3 * self.lambdaHatTheta1Theta3bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(2) * (((self.MTheta1)**(2) + -1 * (self.MTheta3)**(2)))**(-2) * self.yu[fla,flb] + -5/8 * self.lambdaHatTheta1Theta3 * self.lambdaHatTheta1Theta3bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(2) * (((self.MTheta1)**(2) + -1 * (self.MTheta3)**(2)))**(-2) * self.yu[fla,flb] + 15/16 * self.lambdaHatTheta1Theta3 * self.lambdaHatTheta1Theta3bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(4) * (((self.MTheta1)**(2) + -1 * (self.MTheta3)**(2)))**(-3) * self.yu[fla,flb] + -5/8 * self.lambdaHatTheta1Theta3 * self.lambdaHatTheta1Theta3bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(4) * (((self.MTheta1)**(2) + -1 * (self.MTheta3)**(2)))**(-3) * np.log((self.MTheta1)**(2) * (self.mu)**(-2)) * self.yu[fla,flb] + 5/8 * self.lambdaHatTheta1Theta3 * self.lambdaHatTheta1Theta3bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(2) * (((self.MTheta1)**(2) + -1 * (self.MTheta3)**(2)))**(-2) * np.log((self.MTheta1)**(2) * (self.mu)**(-2)) * self.yu[fla,flb] + 5/8 * self.lambdaHatTheta1Theta3 * self.lambdaHatTheta1Theta3bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(2) * (((self.MTheta1)**(2) + -1 * (self.MTheta3)**(2)))**(-2) * np.log((self.MTheta3)**(2) * (self.mu)**(-2)) * self.yu[fla,flb] + 5/8 * self.lambdaHatTheta1Theta3 * self.lambdaHatTheta1Theta3bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(4) * (((self.MTheta1)**(2) + -1 * (self.MTheta3)**(2)))**(-3) * np.log((self.MTheta3)**(2) * (self.mu)**(-2)) * self.yu[fla,flb])

    def alphaOdH(self, fla,flb):
        return (-1/384 * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) * (self.g2)**(4) * self.yd[fla,flb] + -1/384 * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(-2) * (self.g2)**(4) * self.yd[fla,flb] + 5/96 * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) * (self.lambdaHatPrimeTheta1)**(2) * self.yd[fla,flb] + 5/96 * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(-2) * (self.lambdaHatPrimeTheta3)**(2) * self.yd[fla,flb] + 1/18 * self.lambdaHatPrimePrimeTheta1 * self.lambdaHatPrimePrimeTheta1bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) * self.yd[fla,flb] + 3/16 * self.lambdaTheta3 * self.lambdaTheta3bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(-2) * self.yd[fla,flb] + 5/48 * self.lambdaTheta1 * self.lambdaTheta1bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) * self.yd[fla,flb] + -15/16 * self.lambdaHatTheta1Theta3 * self.lambdaHatTheta1Theta3bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(4) * (((self.MTheta1)**(2) + -1 * (self.MTheta3)**(2)))**(-3) * self.yd[fla,flb] + -5/8 * self.lambdaHatTheta1Theta3 * self.lambdaHatTheta1Theta3bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(2) * (((self.MTheta1)**(2) + -1 * (self.MTheta3)**(2)))**(-2) * self.yd[fla,flb] + -5/8 * self.lambdaHatTheta1Theta3 * self.lambdaHatTheta1Theta3bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(2) * (((self.MTheta1)**(2) + -1 * (self.MTheta3)**(2)))**(-2) * self.yd[fla,flb] + 15/16 * self.lambdaHatTheta1Theta3 * self.lambdaHatTheta1Theta3bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(4) * (((self.MTheta1)**(2) + -1 * (self.MTheta3)**(2)))**(-3) * self.yd[fla,flb] + -5/8 * self.lambdaHatTheta1Theta3 * self.lambdaHatTheta1Theta3bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(4) * (((self.MTheta1)**(2) + -1 * (self.MTheta3)**(2)))**(-3) * np.log((self.MTheta1)**(2) * (self.mu)**(-2)) * self.yd[fla,flb] + 5/8 * self.lambdaHatTheta1Theta3 * self.lambdaHatTheta1Theta3bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(2) * (((self.MTheta1)**(2) + -1 * (self.MTheta3)**(2)))**(-2) * np.log((self.MTheta1)**(2) * (self.mu)**(-2)) * self.yd[fla,flb] + 5/8 * self.lambdaHatTheta1Theta3 * self.lambdaHatTheta1Theta3bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(2) * (((self.MTheta1)**(2) + -1 * (self.MTheta3)**(2)))**(-2) * np.log((self.MTheta3)**(2) * (self.mu)**(-2)) * self.yd[fla,flb] + 5/8 * self.lambdaHatTheta1Theta3 * self.lambdaHatTheta1Theta3bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(4) * (((self.MTheta1)**(2) + -1 * (self.MTheta3)**(2)))**(-3) * np.log((self.MTheta3)**(2) * (self.mu)**(-2)) * self.yd[fla,flb])

    def alphaOeH(self, fla,flb):
        return (-1/384 * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) * (self.g2)**(4) * self.yl[fla,flb] + -1/384 * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(-2) * (self.g2)**(4) * self.yl[fla,flb] + 5/96 * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) * (self.lambdaHatPrimeTheta1)**(2) * self.yl[fla,flb] + 5/96 * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(-2) * (self.lambdaHatPrimeTheta3)**(2) * self.yl[fla,flb] + 1/18 * self.lambdaHatPrimePrimeTheta1 * self.lambdaHatPrimePrimeTheta1bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) * self.yl[fla,flb] + 3/16 * self.lambdaTheta3 * self.lambdaTheta3bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(-2) * self.yl[fla,flb] + 5/48 * self.lambdaTheta1 * self.lambdaTheta1bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) * self.yl[fla,flb] + -15/16 * self.lambdaHatTheta1Theta3 * self.lambdaHatTheta1Theta3bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(4) * (((self.MTheta1)**(2) + -1 * (self.MTheta3)**(2)))**(-3) * self.yl[fla,flb] + -5/8 * self.lambdaHatTheta1Theta3 * self.lambdaHatTheta1Theta3bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(2) * (((self.MTheta1)**(2) + -1 * (self.MTheta3)**(2)))**(-2) * self.yl[fla,flb] + -5/8 * self.lambdaHatTheta1Theta3 * self.lambdaHatTheta1Theta3bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(2) * (((self.MTheta1)**(2) + -1 * (self.MTheta3)**(2)))**(-2) * self.yl[fla,flb] + 15/16 * self.lambdaHatTheta1Theta3 * self.lambdaHatTheta1Theta3bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(4) * (((self.MTheta1)**(2) + -1 * (self.MTheta3)**(2)))**(-3) * self.yl[fla,flb] + -5/8 * self.lambdaHatTheta1Theta3 * self.lambdaHatTheta1Theta3bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(4) * (((self.MTheta1)**(2) + -1 * (self.MTheta3)**(2)))**(-3) * np.log((self.MTheta1)**(2) * (self.mu)**(-2)) * self.yl[fla,flb] + 5/8 * self.lambdaHatTheta1Theta3 * self.lambdaHatTheta1Theta3bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(2) * (((self.MTheta1)**(2) + -1 * (self.MTheta3)**(2)))**(-2) * np.log((self.MTheta1)**(2) * (self.mu)**(-2)) * self.yl[fla,flb] + 5/8 * self.lambdaHatTheta1Theta3 * self.lambdaHatTheta1Theta3bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(2) * (((self.MTheta1)**(2) + -1 * (self.MTheta3)**(2)))**(-2) * np.log((self.MTheta3)**(2) * (self.mu)**(-2)) * self.yl[fla,flb] + 5/8 * self.lambdaHatTheta1Theta3 * self.lambdaHatTheta1Theta3bar * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(4) * (((self.MTheta1)**(2) + -1 * (self.MTheta3)**(2)))**(-3) * np.log((self.MTheta3)**(2) * (self.mu)**(-2)) * self.yl[fla,flb])

    def alphaOuGbar(self, ):
        return 0

    def alphaOuWbar(self, ):
        return 0

    def alphaOuBbar(self, ):
        return 0

    def alphaOdGbar(self, ):
        return 0

    def alphaOdWbar(self, ):
        return 0

    def alphaOdBbar(self, ):
        return 0

    def alphaOeWbar(self, ):
        return 0

    def alphaOeBbar(self, ):
        return 0

    def alphaOHq1bar(self, ):
        return 0

    def alphaOHq3bar(self, ):
        return 0

    def alphaOHubar(self, ):
        return 0

    def alphaOHdbar(self, ):
        return 0

    def alphaOHudbar(self, ):
        return 0

    def alphaOHl1bar(self, ):
        return 0

    def alphaOHl3bar(self, ):
        return 0

    def alphaOHebar(self, ):
        return 0

    def alphaOuHbar(self, ):
        return 0

    def alphaOdHbar(self, ):
        return 0

    def alphaOeHbar(self, ):
        return 0

    def alphaOqq1(self, fl93,fl95,fla,flb):
        return (-1/3840 * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(-2) * (self.g1)**(4) * self.kronecker_delta(fl93,fl95) * self.kronecker_delta(fla,flb) + -1/34560 * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) * (self.g1)**(4) * self.kronecker_delta(fl93,fl95) * self.kronecker_delta(fla,flb))

    def alphaOqq3(self, fl93,fl95,fla,flb):
        return (-1/768 * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) * (self.g2)**(4) * self.kronecker_delta(fl93,fl95) * self.kronecker_delta(fla,flb) + -1/768 * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(-2) * (self.g2)**(4) * self.kronecker_delta(fl93,fl95) * self.kronecker_delta(fla,flb))

    def alphaOuu(self, fl93,fl95,fla,flb):
        return (-1/240 * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(-2) * (self.g1)**(4) * self.kronecker_delta(fl93,fl95) * self.kronecker_delta(fla,flb) + -1/2160 * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) * (self.g1)**(4) * self.kronecker_delta(fl93,fl95) * self.kronecker_delta(fla,flb))

    def alphaOdd(self, fl93,fl95,fla,flb):
        return (-1/960 * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(-2) * (self.g1)**(4) * self.kronecker_delta(fl93,fl95) * self.kronecker_delta(fla,flb) + -1/8640 * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) * (self.g1)**(4) * self.kronecker_delta(fl93,fl95) * self.kronecker_delta(fla,flb))

    def alphaOud1(self, fl93,fl95,fla,flb):
        return (1/240 * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(-2) * (self.g1)**(4) * self.kronecker_delta(fl93,fl95) * self.kronecker_delta(fla,flb) + 1/2160 * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) * (self.g1)**(4) * self.kronecker_delta(fl93,fl95) * self.kronecker_delta(fla,flb))

    def alphaOud8(self, ):
        return 0

    def alphaOqu1(self, fla,flb,fl95,fl93):
        return (-1/480 * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(-2) * (self.g1)**(4) * self.kronecker_delta(fl93,fl95) * self.kronecker_delta(fla,flb) + -1/4320 * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) * (self.g1)**(4) * self.kronecker_delta(fl93,fl95) * self.kronecker_delta(fla,flb))

    def alphaOqu8(self, ):
        return 0

    def alphaOqd1(self, fl95,fl93,fla,flb):
        return (1/960 * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(-2) * (self.g1)**(4) * self.kronecker_delta(fl93,fl95) * self.kronecker_delta(fla,flb) + 1/8640 * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) * (self.g1)**(4) * self.kronecker_delta(fl93,fl95) * self.kronecker_delta(fla,flb))

    def alphaOqd8(self, ):
        return 0

    def alphaOquqd1(self, ):
        return 0

    def alphaOquqd8(self, ):
        return 0

    def alphaOll(self, fl93,fl95,fla,flb):
        return (-3/1280 * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(-2) * (self.g1)**(4) * self.kronecker_delta(fl93,fl95) * self.kronecker_delta(fla,flb) + -1/384 * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) * (self.g2)**(4) * self.kronecker_delta(fl93,flb) * self.kronecker_delta(fl95,fla) + -1/384 * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(-2) * (self.g2)**(4) * self.kronecker_delta(fl93,flb) * self.kronecker_delta(fl95,fla) + -1/3840 * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) * (self.g1)**(4) * self.kronecker_delta(fl93,fl95) * self.kronecker_delta(fla,flb) + 1/768 * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) * (self.g2)**(4) * self.kronecker_delta(fl93,fl95) * self.kronecker_delta(fla,flb) + 1/768 * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(-2) * (self.g2)**(4) * self.kronecker_delta(fl93,fl95) * self.kronecker_delta(fla,flb))

    def alphaOee(self, fl93,fl95,fla,flb):
        return (-3/320 * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(-2) * (self.g1)**(4) * self.kronecker_delta(fl93,fl95) * self.kronecker_delta(fla,flb) + -1/960 * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) * (self.g1)**(4) * self.kronecker_delta(fl93,fl95) * self.kronecker_delta(fla,flb))

    def alphaOle(self, fl95,fl93,fla,flb):
        return (-3/320 * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(-2) * (self.g1)**(4) * self.kronecker_delta(fl93,fl95) * self.kronecker_delta(fla,flb) + -1/960 * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) * (self.g1)**(4) * self.kronecker_delta(fl93,fl95) * self.kronecker_delta(fla,flb))

    def alphaOlq1(self, fla,flb,fl93,fl95):
        return (1/640 * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(-2) * (self.g1)**(4) * self.kronecker_delta(fl93,fl95) * self.kronecker_delta(fla,flb) + 1/5760 * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) * (self.g1)**(4) * self.kronecker_delta(fl93,fl95) * self.kronecker_delta(fla,flb))

    def alphaOlq3(self, fla,flb,fl93,fl95):
        return (-1/384 * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) * (self.g2)**(4) * self.kronecker_delta(fl93,fl95) * self.kronecker_delta(fla,flb) + -1/384 * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(-2) * (self.g2)**(4) * self.kronecker_delta(fl93,fl95) * self.kronecker_delta(fla,flb))

    def alphaOeu(self, fl93,fl95,fla,flb):
        return (1/80 * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(-2) * (self.g1)**(4) * self.kronecker_delta(fl93,fl95) * self.kronecker_delta(fla,flb) + 1/720 * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) * (self.g1)**(4) * self.kronecker_delta(fl93,fl95) * self.kronecker_delta(fla,flb))

    def alphaOed(self, fl93,fl95,fla,flb):
        return (-1/160 * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(-2) * (self.g1)**(4) * self.kronecker_delta(fl93,fl95) * self.kronecker_delta(fla,flb) + -1/1440 * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) * (self.g1)**(4) * self.kronecker_delta(fl93,fl95) * self.kronecker_delta(fla,flb))

    def alphaOqe(self, fl95,fl93,fla,flb):
        return (1/320 * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(-2) * (self.g1)**(4) * self.kronecker_delta(fl93,fl95) * self.kronecker_delta(fla,flb) + 1/2880 * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) * (self.g1)**(4) * self.kronecker_delta(fl93,fl95) * self.kronecker_delta(fla,flb))

    def alphaOlu(self, fla,flb,fl95,fl93):
        return (1/160 * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(-2) * (self.g1)**(4) * self.kronecker_delta(fl93,fl95) * self.kronecker_delta(fla,flb) + 1/1440 * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) * (self.g1)**(4) * self.kronecker_delta(fl93,fl95) * self.kronecker_delta(fla,flb))

    def alphaOld(self, fl95,fl93,fla,flb):
        return (-1/320 * self.onelooporder * (np.pi)**(-2) * (self.MTheta3)**(-2) * (self.g1)**(4) * self.kronecker_delta(fl93,fl95) * self.kronecker_delta(fla,flb) + -1/2880 * self.onelooporder * (np.pi)**(-2) * (self.MTheta1)**(-2) * (self.g1)**(4) * self.kronecker_delta(fl93,fl95) * self.kronecker_delta(fla,flb))

    def alphaOledq(self, ):
        return 0

    def alphaOlequ1(self, ):
        return 0

    def alphaOlequ3(self, ):
        return 0

    def alphaOqq1bar(self, ):
        return 0

    def alphaOqq3bar(self, ):
        return 0

    def alphaOuubar(self, ):
        return 0

    def alphaOddbar(self, ):
        return 0

    def alphaOud1bar(self, ):
        return 0

    def alphaOud8bar(self, ):
        return 0

    def alphaOqu1bar(self, ):
        return 0

    def alphaOqu8bar(self, ):
        return 0

    def alphaOqd1bar(self, ):
        return 0

    def alphaOqd8bar(self, ):
        return 0

    def alphaOquqd1bar(self, ):
        return 0

    def alphaOquqd8bar(self, ):
        return 0

    def alphaOllbar(self, ):
        return 0

    def alphaOeebar(self, ):
        return 0

    def alphaOlebar(self, ):
        return 0

    def alphaOlq1bar(self, ):
        return 0

    def alphaOlq3bar(self, ):
        return 0

    def alphaOeubar(self, ):
        return 0

    def alphaOedbar(self, ):
        return 0

    def alphaOqebar(self, ):
        return 0

    def alphaOlubar(self, ):
        return 0

    def alphaOldbar(self, ):
        return 0

    def alphaOledqbar(self, ):
        return 0

    def alphaOlequ1bar(self, ):
        return 0

    def alphaOlequ3bar(self, ):
        return 0

    def alphaWeinberg(self, ):
        return 0

    def alphaWeinbergbar(self, ):
        return 0
