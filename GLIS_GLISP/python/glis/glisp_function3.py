from numpy import sum as npsum
from numpy import vstack, append, array, size, max


class glisp_function3:
    """
    Preference query function.

    Values of f already computed are stored to save computations
    of expensive functions.

    (C) 2019 by A. Bemporad, September 23, 2019

    Note: Modified to account for unknown constraints (feasibility and satisfactory constraints) by M. Zhu, May 31, 2021
    """

    def __init__(self, f, comparetol, g_unkn, s_unkn):
        self.itest = 0
        self.Xtest = []
        self.Ftest = []
        self.Gtest = []
        self.Stest = []
        self.Festest = []
        self.SoftconstTest = []
        self.f = f
        self.comparetol = comparetol
        self.g_unkn = g_unkn
        self.s_unkn = s_unkn

    def clear(self):
        self.itest = 0
        self.Xtest = []
        self.Ftest = []
        self.Gtest = []
        self.Stest = []
        self.Festest = []
        self.SoftconstTest = []
        return

    def eval(self, x, y):

        xfound = False
        yfound = False

        itest = self.itest
        if itest > 0:
            Xtest = self.Xtest
            Ftest = self.Ftest
            Gtest = self.Gtest
            Stest = self.Stest
            Festest = self.Festest
            SoftconstTest = self.SoftconstTest
        else:
            fx = self.f(x)
            gx_unkn = self.g_unkn(x)
            sx_unkn = self.s_unkn(x)
            if gx_unkn < self.comparetol:
                fesx = 1
            else:
                fesx = 0
            if sx_unkn < self.comparetol:
                softx = 1
            else:
                softx =0

            itest = 1
            Xtest = array([x])
            Ftest = array([fx])
            Gtest = array([gx_unkn])
            Stest = array([sx_unkn])
            Festest = array([fesx])
            SoftconstTest = array([softx])
            xfound = True

        for i in range(itest):
            if not (xfound) and npsum(abs(Xtest[i, :] - x)) <= 1e-10:
                xfound = True
                fx = Ftest[i]
                fesx = Festest[i]
                softx = SoftconstTest[i]
                gx_unkn = Gtest[i]
                sx_unkn= Stest[i]

            if not (yfound) and npsum(abs(Xtest[i, :] - y)) <= 1e-10:
                yfound = True
                fy = Ftest[i]
                fesy = Festest[i]
                softy = SoftconstTest[i]
                gy_unkn = Gtest[i]
                sy_unkn = Stest[i]

        if not (xfound):
            fx = self.f(x)
            gx_unkn = self.g_unkn(x)
            sx_unkn = self.s_unkn(x)
            if gx_unkn < self.comparetol:
                fesx = 1
            else:
                fesx = 0
            if sx_unkn < self.comparetol:
                softx = 1
            else:
                softx = 0

            Xtest = vstack((Xtest, x))
            Ftest = append(Ftest, fx)
            Gtest = append(Gtest,gx_unkn)
            Stest = append(Stest,sx_unkn)
            Festest = append(Festest, fesx)
            SoftconstTest = append(SoftconstTest,softx)
            itest = itest + 1

        if not (yfound):
            fy = self.f(y)
            gy_unkn = self.g_unkn(y)
            sy_unkn = self.s_unkn(y)
            if gy_unkn < self.comparetol:
                fesy = 1
            else:
                fesy = 0
            if sy_unkn < self.comparetol:
                softy =1
            else:
                softy = 0
            Xtest = vstack((Xtest, y))
            Ftest = append(Ftest, fy)
            Gtest = append(Gtest,gy_unkn)
            Stest = append(Stest,sy_unkn)
            Festest = append(Festest, fesy)
            SoftconstTest = append(SoftconstTest,softy)
            itest = itest + 1

        # Make comparison
        if fx < fy - self.comparetol:
            if (fesx+softx ==2) or (fesx+softx ==0 and fesy+softy ==0):
                out =-1
            else:
                out = 1
        elif fx > fy - self.comparetol:
            if (fesy+softy ==2) or (fesx+softx ==0 and fesy+softy ==0):
                out = 1
            else:
                out =-1
        else:
            out =0


        self.Xtest = Xtest
        self.Ftest = Ftest
        self.Gtest = Gtest
        self.Stest = Stest
        self.itest = itest
        self.Festest = Festest
        self.SoftconstTest = SoftconstTest

        return out, fesx, fesy, softx, softy

    def value(self, x):
        # Compute function value & the feasibility, from available ones if available
        #
        # (C) 2019 A. Bemporad, September 22, 2019
        # modified May, 2021, MZ

        j = 0
        while j < self.itest:
            if npsum(abs(self.Xtest[j, :] - x)) <= 1e-10:
                val = self.Ftest[j]
                gx_unkn = self.Gtest[j]
                fes = self.Festest[j]
                return val, fes
            j = j + 1

        # Value does not exist, compute it
        val = self.f(x)
        gx_unkn = self.g_unkn(x)
        sx_unkn = self.s_unkn(x)
        if gx_unkn < self.comparetol:
            fes = 1
        else:
            fes = 0
        if sx_unkn < self.comparetol:
            soft =1
        else:
            soft =0
        self.Xtest = vstack((self.Xtest, x))
        self.Ftest = append(self.Ftest, val)
        self.Gtest = append(self.Gtest,gx_unkn)
        self.Stest = append(self.Stest,sx_unkn)
        self.Festest = append(self.Festest, fes)
        self.SoftconstTest = append(self.SoftconstTest,soft)
        self.itest = self.itest + 1
        return val, fes, soft



