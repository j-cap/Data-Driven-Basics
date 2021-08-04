def solve(prob):
    """
    Solve (GL)obal optimization problems using (I)nverse distance
    weighting and radial basis function (S)urrogates.

    (C) 2019 A. Bemporad

    sol = glis.solve(prob) solves the global optimization problem

    min  f(x)
    s.t. lb <= x <=ub, A*x <=b, g(x)<=0

    using the global optimization algorithm described in [1]. The approach is
    particularly useful when f(x) is time-consuming to evaluate, as it
    attempts at minimizing the number of function evaluations.

    The default problem structure is

    prob = glis.default(nvars)

    where nvars = dimension of optimization vector x. See function glis_default
    for a description of all available options.

    The output argument 'out' is a structure reporting the following information:

    out["X"]:    trace of all samples x at which f(x) has been evaluated
    out["F"]:    trace of all function evaluations f(x)
    out["W"]:    final set of weights (only meaningful for RBFs)
    out["M"]:    RBF matrix (only meaningful for RBFs)
    out["xopt"]: best sample found during search
    out["fopt"]: best value found during search, fopt=f(xopt)
	
	out["time_iter"], out["time_opt_acquisition"], out["time_fit_surrogate"], out["time_f_eval"]
	store timing recorded during the execution of the algorithm.

    Required Python packages:
        pyDOE:   https://pythonhosted.org/pyDOE/
        nlopt:   https://nlopt.readthedocs.io (required only if DIRECT solver is used)
        pyswarm: https://pythonhosted.org/pyswarm/ (required only if PSO solver is used)

    [1] A. Bemporad, "Global optimization via inverse weighting and radial basis functions,"
        Computational Optimization and Applications, vol. 77, pp. 571–595.


    # (C-GLIS)
    # Note: Add features to handle unknown constraints (M. Zhu, June 08, 2021)
    #       Known constraints will be handled via penalty functions
    #       For unknown constraints, here we assume that after performing the experiment, we can access the feasibility & satisfactory labels
    #
    # Following are the new parameters introduced in C-GLISp
    # opts["isUnknownFeasibilityConstrained"]: if true, unknown feasibility constraints are involved
    # opts["isUnknownSatisfactionConstrained"]: if true, unknown satisfaction constraints are involed
    # delta_E: delta for te pure IDW exploration term, \delta_E in the paper
    # delta_G_default: delta for feasibility constraints, \delta_{G,default} in the paper
    # delta_S_default: delta for satisfaction constraints, \delta_{S,default} in the paper
    # Feasibility_unkn: feasibility labels for unknown feasibility constraints
    # SatConst_unkn: satisfaction labels for unknown satisfactory constraints
    """

    import glis.glis_init as glis_init
    from pyswarm import pso  # https://pythonhosted.org/pyswarm/

    from numpy.linalg import svd

    from numpy import zeros, ones, diag, inf, argmin
    from numpy import where, maximum, exp, vstack
    from numpy import sum as npsum
    from numpy import array as nparray
    from math import sqrt, atan, pi
    import contextlib
    import io
    import time

    def get_rbf_weights(M, F, NX, svdtol):
        # Solve M*W = F using SVD

        U, dS, V = svd(M[0:NX, 0:NX])
        ii = where(dS >= svdtol)
        ns = max(ii[0]) + 1
        W = (V[0:ns, ].T).dot(diag(1 / dS[0:ns].flatten('C')).dot((U[:, 0:ns].T).dot(F[0:NX])))

        return W

    def get_delta_adpt(X,constraint_set,delta_const_default):
        ind = constraint_set.shape[0]
        sqr_error_feas = zeros((ind,1))
        for i in range(0,ind):
            xx = X[i,:]
            Xi = vstack((X[0:i, :], X[i + 1:ind, :]))
            constraint_set_i = vstack((constraint_set[0:i,],constraint_set[i+1:ind,]))
            Feas_xx = constraint_set[i]
            d = npsum((Xi - xx) ** 2, axis=-1)
            w = npsum(-d)/d
            sw = sum(w)
            ghat = npsum(constraint_set_i.T * w) / sw
            sqr_error_feas[i] = (ghat-Feas_xx)**2

        std_feas = (sum(sqr_error_feas)/(ind-1))**(1/2)
        delta_adpt = (1-std_feas) *delta_const_default

        return delta_adpt

    def facquisition(xx, X, F, N, alpha, delta_E, dF, W, rbf, useRBF, isUnknownFeasibilityConstrained,isUnknownSatisfactionConstrained,Feasibility_unkn,SatConst_unkn,delta_G,delta_S,iw_ibest,maxevals):
        # Acquisition function to minimize to get next sample

        d = npsum((X[0:N, ] - xx) ** 2, axis=-1)

        ii = where(d < 1e-12)
        if ii[0].size > 0:
            fhat = F[ii[0]][0]
            dhat = 0
            if isUnknownFeasibilityConstrained:
                Ghat = Feasibility_unkn[ii]
            else:
                Ghat = 1
            if isUnknownSatisfactionConstrained:
                Shat = SatConst_unkn[ii]
            else:
                Shat = 1
        else:
            w = exp(-d) / d
            sw = sum(w)

            if useRBF:
                v = rbf(X[0:N,:],xx)
                fhat = v.ravel().dot(W.ravel())
            else:
                fhat = npsum(F[0:N, ] * w) / sw

            if maxevals <= 30:
  				# for comparision, used in the original GLIS and when N_max <= 30 in C-GLIS
                dhat = delta_E * atan(1 / sum(1 / d)) * 2 / pi * dF + alpha * sqrt(sum(w * (F[0:N, ] - fhat).flatten("c") ** 2) / sw)  
            else:
                dhat = delta_E *((1 - N / maxevals) * atan((1 / sum(1. / d)) / iw_ibest) + N / maxevals * atan(1 / sum(1. / d)))* 2 / pi * dF + alpha * sqrt(sum(w * (F[0:N, ] - fhat).flatten("c") ** 2) / sw)

            # to account for the unknown constraints
            if isUnknownFeasibilityConstrained:
                Ghat = npsum(Feasibility_unkn[0:N].T * w) / sw
            else:
                Ghat = 1

            if isUnknownSatisfactionConstrained:
                Shat = npsum(SatConst_unkn[0:N].T * w) / sw
            else:
                Shat = 1

        f = fhat - dhat + (delta_G*(1-Ghat)+delta_S*(1-Shat))*dF

        return f

    (f, lb, ub, nvar, Aineq, bineq, g, isLinConstrained, isNLConstrained,
     X, F, z, nsamp, maxevals, epsDeltaF, alpha, delta, rhoC, display, svdtol,
     dd, d0, useRBF, rbf, M, scalevars, globoptsol, DIRECTopt,
     PSOiters, PSOswarmsize) = glis_init.init(prob)


    time_iter = []
    time_f_eval = []
    time_opt_acquisition = []
    time_fit_surrogate = []

    for i in range(nsamp):
        time_fun_eval_start = time.perf_counter()
        F[i] = f(X[i,].T)
        time_fun_eval_i = time.perf_counter() - time_fun_eval_start
        time_iter.append(time_fun_eval_i)
        time_f_eval.append(time_fun_eval_i)
        time_opt_acquisition.append(0.0)
        time_fit_surrogate.append(0.0)

    isUnknownFeasibilityConstrained = prob["isUnknownFeasibilityConstrained"]
    isUnknownSatisfactionConstrained = prob['isUnknownSatisfactionConstrained']
    g0 = prob["g_unkn_fun"]
    s0 = prob["s_unkn_fun"]
    if scalevars:
        g_unkn = lambda x: g0(x * dd + d0)
        s_unkn = lambda x: s0(x * dd + d0)

    delta_E = delta
    delta_G_default = delta
    delta_S_default = delta/2
    Feasibility_unkn = [] # feasibility labels
    SatConst_unkn = [] # satisfactory constraint label
    isfeas_seq = ones((maxevals,1)).astype(int) # keep track the feasibility of the decision variables(including both known and unknown constraints)
    ibestseq = ones((maxevals,1)).astype(int) # keep track of the ibest throughout

    # Initial sampling phase
    if isUnknownFeasibilityConstrained:
        Feasibility_unkn = zeros((nsamp, 1))
        for i in range(1, nsamp):
            Feasibility_unkn[i] = g_unkn(X[i,:]) < 1e-6
        delta_G = get_delta_adpt(X, Feasibility_unkn, delta_G_default)
    else:
        Feasibility_unkn = ones((nsamp, 1))
        delta_G = 0

    if isUnknownSatisfactionConstrained:
        SatConst_unkn = zeros((nsamp, 1))
        for i in range(1, nsamp):
            SatConst_unkn[i] = s_unkn(X[i,:]) < 1e-6
        delta_S = get_delta_adpt(X, SatConst_unkn, delta_S_default)
    else:
        SatConst_unkn = ones((nsamp, 1))
        delta_S = 0

    if useRBF:
        W = get_rbf_weights(M, F, nsamp, svdtol)
    else:
        W = []

    fbest = inf
    ibest = argmin(F)
    zbest = X[ibest,:]
    # zbest = zeros((nsamp, 1))
    for i in range(nsamp):
        isfeas = True
        if isLinConstrained:
            isfeas = isfeas and all(Aineq.dot(X[i,].T) <= bineq.flatten("c"))
        if isNLConstrained:
            isfeas = isfeas and all(g(X[i,]) <= 0)
        if isUnknownFeasibilityConstrained:
            isfeas = isfeas and Feasibility_unkn[i]> 0
        if isUnknownSatisfactionConstrained:
            isfeas = isfeas and SatConst_unkn[i] > 0
        if isfeas and fbest > F[i]:
            fbest = F[i]
            zbest = X[i,]
            ibest = i

        ibestseq[i] = ibest
        isfeas_seq[i] = isfeas

    Fmax = max(F[0:nsamp])
    Fmin = min(F[0:nsamp])

    N = nsamp

    while N < maxevals:

        time_iter_start = time.perf_counter()

        dF = Fmax - Fmin

        if isLinConstrained or isNLConstrained:
            penalty = rhoC * dF
        if isLinConstrained and isNLConstrained:
            constrpenalty = lambda x: (penalty * (sum(maximum((Aineq.dot(x) - bineq).flatten("c"), 0) ** 2)
                                                  + sum(maximum(g(x), 0) ** 2)))
        elif isLinConstrained and not isNLConstrained:
            constrpenalty = lambda x: penalty * (sum(maximum((Aineq.dot(x) - bineq).flatten("c"), 0) ** 2))
        elif not isLinConstrained and isNLConstrained:
            constrpenalty = lambda x: penalty * sum(maximum(g(x), 0) ** 2)
        else:
            constrpenalty = lambda x: 0

        d_ibest = npsum((vstack((X[0:ibest,:],X[ibest+1:,:]))-X[ibest,:])**2,axis=-1)
        ii = where(d_ibest<1e-12)
        if ii[0].size > 0:
            iw_ibest = 0
        else:
            iw_ibest = 1.0/sum(1.0/d_ibest)

        acquisition = lambda x: (facquisition(x, X, F, N, alpha, delta, dF, W, rbf, useRBF, isUnknownFeasibilityConstrained,isUnknownSatisfactionConstrained,Feasibility_unkn,SatConst_unkn,delta_G,delta_S,iw_ibest,maxevals)
                                 + constrpenalty(x))

        time_opt_acq_start = time.perf_counter()
        if globoptsol == "pswarm":
            # pso(func, lb, ub, ieqcons=[], f_ieqcons=None, args=(), kwargs={},
            #    swarmsize=100, omega=0.5, phip=0.5, phig=0.5, maxiter=100, minstep=1e-8,
            #    minfunc=1e-8, debug=False)
            with contextlib.redirect_stdout(io.StringIO()):
                z, cost = pso(acquisition, lb, ub, swarmsize=PSOswarmsize,
                              minfunc=dF * 1e-8, maxiter=PSOiters)

        elif globoptsol == "direct":
            DIRECTopt.set_min_objective(lambda x, grad: acquisition(x)[0])
            z = DIRECTopt.optimize(z.flatten("c"))
        time_opt_acquisition.append(time.perf_counter() - time_opt_acq_start)

        time_fun_eval_start = time.perf_counter()
        fz = f(z)  # function evaluation
        time_f_eval.append(time.perf_counter() - time_fun_eval_start)

        N = N + 1

        X[N - 1,] = z.T
        F[N - 1] = fz

        Fmax = max(Fmax, fz)
        Fmin = min(Fmin, fz)

        if isUnknownFeasibilityConstrained:
            fesN = g_unkn(z) < 1e-6
            Feasibility_unkn = vstack((Feasibility_unkn, fesN))
            delta_G = get_delta_adpt(X, Feasibility_unkn, delta_G_default)
        else:
            delta_G = 0

        if isUnknownSatisfactionConstrained:
            satconstN = s_unkn(z) < 1e-6
            SatConst_unkn = vstack((SatConst_unkn, satconstN))
            delta_S = get_delta_adpt(X, SatConst_unkn, delta_S_default)
        else:
            delta_S = 0

        isfeas = True
        if isLinConstrained:
            isfeas = isfeas and all(Aineq.dot(z) <= bineq.flatten("c"))
        if isNLConstrained:
            isfeas = isfeas and all(g(z) <= 0)
        if isUnknownFeasibilityConstrained:
            isfeas = isfeas and Feasibility_unkn[N-1]> 0
        if isUnknownSatisfactionConstrained:
            isfeas = isfeas and SatConst_unkn[N-1] > 0
        if isfeas and fbest > fz:
            fbest = fz.copy()
            zbest = z.copy()
            ibest = N-1

        ibestseq[N-1] = ibest
        isfeas_seq[N-1] = isfeas


        time_fit_surrogate_start = time.perf_counter()
        if useRBF:
            # Just update last row and column of M
            for h in range(N):
                mij = rbf(X[h,], X[N - 1,])
                M[h, N - 1] = mij
                M[N - 1, h] = mij

            W = get_rbf_weights(M, F, N, svdtol)
        
        time_fit_surrogate.append(time.perf_counter() - time_fit_surrogate_start)

        if display > 0:

            print("N = %4d, cost = %7.4f, best = %7.4f" % (N, fz, fbest))

            string = ""
            for j in range(nvar):
                aux = zbest[j]
                if scalevars:
                    aux = aux * dd[j] + d0[j]

                string = string + " x" + str(j + 1) + " = " + ('%7.4f' % aux)
            print(string)

        time_iter.append(time.perf_counter() - time_iter_start)

    # end while
    if isfeas_seq[ibest] == 0:  # for the case where no feasible optimizer is identified
        fbest = min(F)
        ibest = argmin(F)
        zbest=X[ibest,:]
    xopt = zbest.copy()
    if not isUnknownFeasibilityConstrained:
        Feasibility_unkn = ones((maxevals,1))
    if not isUnknownSatisfactionConstrained and not isUnknownFeasibilityConstrained:
        SatConst_unkn = ones((maxevals,1))
    elif not isUnknownSatisfactionConstrained and isUnknownFeasibilityConstrained:
        SatConst_unkn = Feasibility_unkn
    fes_opt_unkn = Feasibility_unkn[ibest]
    satConst_opt_unkn = SatConst_unkn[ibest]
    feas_opt_comb = isfeas_seq[ibest]

    if scalevars:
        # Scale variables back
        xopt = xopt * dd + d0
        X = X * (ones((N, 1)) * dd) + ones((N, 1)) * d0

    fopt = fbest.copy()

    if not useRBF:
        W = []

    out = {"xopt": xopt,
           "fopt": fopt,
           "X": X,
           "F": F,
           "W": W,
           "Feasibility_unkn": Feasibility_unkn,
           "fes_opt_unkn": fes_opt_unkn,
           "SatConst_unkn": SatConst_unkn,
           "satConst_opt_unkn": satConst_opt_unkn,
           "isfeas_seq": isfeas_seq,
           "feas_opt_comb": feas_opt_comb,
           "ibest": ibest,
           "ibestseq": ibestseq,
           "time_iter": nparray(time_iter),
           "time_opt_acquisition": nparray(time_opt_acquisition),
           "time_fit_surrogate": nparray(time_fit_surrogate),
           "time_f_eval": nparray(time_f_eval)}

    return out


def default(nvars):
    """ Generate default problem structure for IDW-RBF Global Optimization.

     problem=glis.default(n) generate a default problem structure for a
     an optimization with n variables.

     (C) 2019 by A. Bemporad.
    """

    import glis.glis_default as glis_default
    problem = glis_default.set(nvars)

    return problem
