from numpy import sum as npsum
from numpy import vstack, append, array, size, max


class glisp_function2:
    """ 
    Preference query function.

    Values of f already computed are stored to save computations
    of expensive functions.

    (C) 2019 by A. Bemporad, September 23, 2019

    Note: Modified to account for unknown constraints (feasibility constraints only) by M. Zhu, May 31, 2021
    """

    def __init__(self, f, comparetol, g_unkn):
        self.itest = 0
        self.Xtest = []
        self.Ftest = []
        self.Gtest = []
        self.Festest = []
        self.f = f
        self.comparetol = comparetol
        self.g_unkn = g_unkn

    def clear(self):
        self.itest = 0
        self.Xtest = []
        self.Ftest = []
        self.Gtest = []
        self.Festest = []
        return

    def eval(self, x, y):

        xfound = False
        yfound = False

        itest = self.itest
        if itest > 0:
            Xtest = self.Xtest
            Ftest = self.Ftest
            Gtest = self.Gtest
            Festest = self.Festest
        else:
            fx = self.f(x)
            gx_unkn = self.g_unkn(x)
            if size(gx_unkn) > 1:
                if max(gx_unkn) < self.comparetol:
                    fesx = 1
                else:
                    fesx = 0
            else:
                if gx_unkn < self.comparetol:
                    fesx = 1
                else:
                    fesx = 0

            itest = 1
            Xtest = array([x])
            Ftest = array([fx])
            Gtest = array([gx_unkn])
            Festest = array([fesx])
            xfound = True

        for i in range(itest):
            if not (xfound) and npsum(abs(Xtest[i, :] - x)) <= 1e-10:
                xfound = True
                fx = Ftest[i]
                fesx = Festest[i]
                gx_unkn = Gtest[i]

            if not (yfound) and npsum(abs(Xtest[i, :] - y)) <= 1e-10:
                yfound = True
                fy = Ftest[i]
                fesy = Festest[i]
                gy_unkn = Gtest[i]

        if not (xfound):
            fx = self.f(x)
            gx_unkn = self.g_unkn(x)
            if size(gx_unkn) > 1:
                if max(gx_unkn) < self.comparetol:
                    fesx = 1
                else:
                    fesx = 0
            else:
                if gx_unkn < self.comparetol:
                    fesx = 1
                else:
                    fesx = 0

            Xtest = vstack((Xtest, x))
            Ftest = append(Ftest, fx)
            Gtest = append(Gtest,gx_unkn)
            Festest = append(Festest, fesx)
            itest = itest + 1

        if not (yfound):
            fy = self.f(y)
            gy_unkn = self.g_unkn(y)
            if size(gy_unkn) > 1:
                if max(gy_unkn) < self.comparetol:
                    fesy = 1
                else:
                    fesy = 0
            else:
                if gy_unkn < self.comparetol:
                    fesy = 1
                else:
                    fesy = 0
            Xtest = vstack((Xtest, y))
            Ftest = append(Ftest, fy)
            Gtest = append(Gtest,gy_unkn)
            Festest = append(Festest, fesy)
            itest = itest + 1

        # Make comparison
        if (fesx ==1 and fesy ==1):
            if fx < fy - self.comparetol:
                out = -1
            elif fx > fy + self.comparetol:
                out = 1
            else:
                out = 0
        elif (fesx ==1 and fesy ==0):
            out = -1
        elif (fesx ==0 and fesy ==1):
            out = 1
        else:
            if gx_unkn < gy_unkn- self.comparetol:
                out = -1
            elif gx_unkn > gy_unkn- self.comparetol:
                out = 1
            else:
                out = 0



        self.Xtest = Xtest
        self.Ftest = Ftest
        self.Gtest = Gtest
        self.itest = itest
        self.Festest = Festest

        return out, fesx, fesy

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
        if size(gx_unkn) > 1:
            if max(gx_unkn) < self.comparetol:
                fes = 1
            else:
                fes = 0
        else:
            if gx_unkn < self.comparetol:
                fes = 1
            else:
                fes = 0
        self.Xtest = vstack((self.Xtest, x))
        self.Ftest = append(self.Ftest, val)
        self.Gtest = append(self.Gtest,gx_unkn)
        self.Festest = append(self.Festest, fes)
        self.itest = self.itest + 1
        return val, fes



