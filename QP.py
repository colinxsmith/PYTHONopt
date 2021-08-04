#!/bin/env python
# Colin July 2006
# Interior Point Methods for linear and quadratic programming
# Homogenous method does not work for QP
import sys
sys.path.insert(0, 'c:\\Users\colin\safeqp64')
from Optn import *
def square(x): return x*x
rootp5 = pow(2, -.5)
def gfunc(x): return min(0.5, square(1.-x))*(1.-x)
def norm(x): return pow(dot(x, x), .5)
def differ(a, b): return float(abs(a-b))/a
eps = epsget()
sqeps = int(1e-12/eps)*eps
def inverseDv(n, D, V):  # Not quicker!!!!! (Unless only one V)
    """Inverse of D + sum(vivi') where D is diagonal and vi are orthogonal"""
    DI = [0]*(n*(n+1)/2)
    for i in range(n):
        DI[(i+3)*i/2] = 1./D[i]
    for l in range(n):
        VD = []
        Sym_mult(n, DI, V[l*n:(l+1)*n], VD)
        bot = pow(1+dot(V[l*n:(l+1)*n], VD), .5)
        VD = [i/bot for i in VD]
        ij = 0
        for i in range(n):
            for j in range(i+1):
                DI[ij] -= VD[i]*VD[j]
                ij += 1
    return DI


def Q2V(n, Q):
    """Decompose Q = sum(vivi') where vi mutually orthogonal"""
    (e, V) = eigen(n, Q)
    for i in range(n):
        p = pow(e[i], .5)
        for j in range(n):
            V[j+n*i] *= p
    return V
class LP:
    def __init__(self, n, m, c, A, b):
        for i in 'n m c A b'.split():
            setattr(self, i, eval(i))
        self.x = [1]*n
        self.s = [1]*n
        self.sign = [1]*n
        self.y = [0]*m
        self.tau = 1
        self.kappa = 1
        self.homo = 0
        self.tridiag = 0

    def presid(self):
        Ax = [0]*self.m
        dmxtmultv(self.n, self.m, self.A, self.x, Ax)
        self.rp = [self.b[i]*self.tau - Ax[i] for i in range(self.m)]

    def dresid(self):
        Ay = [1]*self.n
        dmxtmulv(self.m, self.n, self.A, self.y, Ay)
        #Ay=[dot(self.y,self.A[i*self.m:(i+1)*self.m]) for i in range(self.n)]
        self.rd = [self.c[i]*self.tau - self.s[i] - Ay[i]
                   for i in range(self.n)]

    def calcmu(self):
        mu = [0]
        dmxtmulv(self.n, 1, [self.x[i] for i in range(self.n)], self.s, mu)
        if not self.homo:
            print('Calculated mu', mu[0]/self.n)
            self.mu = mu[0]/self.n
        else:
            print('Calculated mu', (mu[0]+self.tau*self.kappa)/(self.n+1))
            self.mu = (mu[0]+self.tau*self.kappa)/(self.n+1)

    def muresid(self):
        self.rmu = [self.mu - self.x[i]*self.s[i] for i in range(self.n)]
        if self.homo:
            self.hrmu = self.mu-self.kappa*self.tau
            self.rkxy = self.kappa+dot(self.c, self.x)-dot(self.b, self.y)

    def getdelta(self):
        self.getv()
        self.delta = 0.5 * \
            pow(sum([square(self.v[i]-1./self.v[i])
                for i in range(self.n)]), .5)

    def solvePC(self, gamma, corrector=0, dx0=[], ds0=[], dtau0=0, dkappa0=0):
        g1 = 1-gamma
        M = []
        if not self.tridiag:
            order = [0]*self.m
        else:
            Q = [0]*(self.m*self.m)
            DD = [0]*self.m
            EE = [0]*(self.m-1)
        if not corrector:
            for i in range(self.m):  # A(X/S)A'
                for j in range(i+1):
                    ss = 0
                    for k in range(self.n):
                        ss += self.A[k*self.m+i] * \
                            self.A[k*self.m+j]*self.x[k]/self.s[k]
                    M.append(ss)
            back = 0
            if not self.tridiag:
                if self.m != 1:
                    back = dsptrf('U', self.m, M, order)  # factorise
            else:
                if self.m != 1:
                    back = trifact(self.m, M, Q, DD, EE)  # factorise
            if back != 0:
                raise 'Equations have gone singular %d' % back
            if not self.tridiag:
                self.M = M
                self.order = order
            else:
                self.Q = Q
                self.DD = DD
                self.EE = EE
        else:
            if not self.tridiag:
                M = self.M
                order = self.order
            else:
                Q = self.Q
                DD = self.DD
                EE = self.EE
        if not corrector:
            w1 = [(self.x[i]*self.rd[i]*g1 - (self.rmu[i]-g1*self.mu))/self.s[i]
                  for i in range(self.n)]
        else:
            w1 = [(self.x[i]*self.rd[i]*g1 - (self.rmu[i]-g1*self.mu -
                   dx0[i]*ds0[i]))/self.s[i] for i in range(self.n)]
        dy = [0]*self.m
        dmxtmultv(self.n, self.m, self.A, w1, dy)
        dy = [dy[i]+self.rp[i]*g1 for i in range(self.m)]  # rhs
        if not self.tridiag:
            if self.m == 1:
                dy[0] /= M[0]
            else:
                dsptrs('U', self.m, 1, M, order, dy, self.m)  # solve for dy
        else:
            if self.m == 1:
                dy[0] /= M[0]
            else:
                trisolve(self.m, Q, DD, EE, dy)  # solve for dy
        if not self.homo:
            ds = [0]*self.n
            dmxtmulv(self.m, self.n, self.A, dy, ds)
            ds = [self.rd[i]*g1-ds[i] for i in range(self.n)]
            if not corrector:
                dx = [(self.rmu[i]-g1*self.mu-ds[i]*self.x[i])/self.s[i]
                      for i in range(self.n)]
            else:
                dx = [(self.rmu[i]-g1*self.mu-ds[i]*self.x[i]-dx0[i]
                       * ds0[i])/self.s[i] for i in range(self.n)]
            return(dx, ds, dy)
        else:
            db = [0]*self.m
            cx = [self.c[i]*self.x[i]/self.s[i] for i in range(self.n)]
            dmxtmultv(self.n, self.m, self.A, cx, db)
            db = [cx[i]+self.b[i] for i in range(self.m)]
            if not self.tridiag:
                if self.m == 1:
                    db[0] /= M[0]
                else:
                    dsptrs('U', self.m, 1, M, order, db, self.m)
            else:
                if self.m == 1:
                    db[0] /= M[0]
                else:
                    trisolve(self.m, Q, DD, EE, db)
            dx = [0]*self.n
            dmxtmulv(self.m, self.n, self.A, dy, dx)
            dx = [self.x[i]/self.s[i]*dx[i]-w1[i] for i in range(self.n)]
            dc = [0]*self.n
            dmxtmulv(self.m, self.n, self.A, db, dc)
            dc = [self.x[i]/self.s[i]*dc[i]-cx[i] for i in range(self.n)]
            cdx = dot(self.c, dx)
            bdy = dot(self.b, dy)
            cdc = dot(self.c, dc)
            bdb = dot(self.b, db)
            # -(cdx + dtaucdc) + bdy + dtaubdb = kappa + cx-by +(mu-taukappa - kappadtau)/tau
            if not corrector:
                dtau = (cdx-bdy+self.rkxy*g1+(self.hrmu-g1*self.mu) /
                        self.tau)/(bdb-cdc+self.kappa/self.tau)
            else:
                dtau = (cdx-bdy+self.rkxy*g1+(self.hrmu-g1*self.mu -
                        dtau0*dkappa0)/self.tau)/(bdb-cdc+self.kappa/self.tau)
            dx = [dx[i]+dtau*dc[i] for i in range(self.n)]
            dy = [dy[i]+dtau*db[i] for i in range(self.m)]
            if not corrector:
                ds = [(self.rmu[i]-g1*self.mu-dx[i]*self.s[i])/self.x[i]
                      for i in range(self.n)]
            else:
                ds = [(self.rmu[i]-g1*self.mu-dx[i]*self.s[i]-dx0[i]
                       * ds0[i])/self.x[i] for i in range(self.n)]
            if not corrector:
                dkappa = (self.hrmu-g1*self.mu-self.kappa*dtau)/self.tau
            else:
                dkappa = (self.hrmu-g1*self.mu-self.kappa *
                          dtau-dtau0*dkappa0)/self.tau
            return(dx, ds, dy, dtau, dkappa)

    def maxstep(self, dx, ds, dtau=0, dkappa=0):
        ddx = 1
        dds = 1
        dd = 1
        for i in range(self.n):
            if dx[i]*self.sign[i] < 0:
                ddx = min(ddx, -self.x[i]/dx[i])
            if ds[i]*self.sign[i] < 0:
                dds = min(dds, -self.s[i]/ds[i])
        if self.homo:
            if dtau < 0:
                dd = min(dd, -self.tau/dtau)
            if dkappa < 0:
                dd = min(dd, -self.kappa/dkappa)
        if dd != 1:
            print('max step', dd)
        if ddx != 1:
            print('max step x', ddx)
        if dds != 1:
            print('max step s', dds)
        return (dd, dds, ddx)

    def update(self, dx, dy, ds, step=1, dtau=0, dkappa=0):
        for i in range(self.n):
            self.x[i] += (dx[i]*step)
            self.s[i] += (ds[i]*step)
            if self.x[i]*self.sign[i] < 0:
                print('bad x', i, self.x[i], dx[i])
            if self.s[i]*self.sign[i] < 0:
                print('bad s', i, self.s[i], ds[i])
        for i in range(self.m):
            self.y[i] += (dy[i])
        if self.homo:
            self.tau += dtau*step
            self.kappa += dkappa*step

    def update_sep(self, dx, dy, ds, stepx=1, stepy=1, steps=1, dtau=0, dkappa=0, steph=1):
        for i in range(self.n):
            self.x[i] += (dx[i]*stepx)
            self.s[i] += (ds[i]*steps)
            if self.x[i]*self.sign[i] < 0:
                print('bad x', i, self.x[i], dx[i])
            if self.s[i]*self.sign[i] < 0:
                print('bad s', i, self.s[i], ds[i])
        for i in range(self.m):
            self.y[i] += (dy[i]*stepy)
        if self.homo:
            self.tau += dtau*steph
            self.kappa += dkappa*steph

    def results(self):
        print('tau', self.tau)
        print('kappa', self.kappa)
        """print '%20s %20s %20s %20s'%('x','s','dual residual','central path resdiual')
        for i in range(self.n):
            print '%20.8e %20.8e %20.8e %20.8e'%(self.x[i],self.s[i],self.rd[i],self.rmu[i])
        print '%20s %20s'%('y','primal residual')
        for i in range(self.m):
            print '%20.8e %20.8e'%(self.y[i],self.rp[i])"""
        print('Primal\t%20.8f' % (self.primal()/self.tau))
        print('Dual\t%20.8f' % (self.dual()/self.tau))
        print('Complementarity', self.compl(), self.tau*self.kappa)

    def primal(self):
        return sum([self.c[i]*self.x[i] for i in range(self.n)])

    def dual(self):
        return sum([self.b[i]*self.y[i] for i in range(self.m)])

    def compl(self):
        return sum([self.x[i]*self.s[i] for i in range(self.n)])
class QP(LP):
    def __init__(self, n, m, c, A, b, H=[]):
        for i in 'n m c A b'.split():
            setattr(self, i, eval(i))
        self.x = [1]*n
        self.s = [1]*n
        self.sign = [1]*n
        self.y = [0]*m
        self.tau = 1
        self.kappa = 1
        self.homo = 0
        self.tridiag = 0
        self.H = H
        self.Hx = [0]*self.n

    def Hmul(self):
        if self.H != []:
            Sym_mult(self.n, self.H, self.x, self.Hx)
            self.cmod = [self.c[i]+self.Hx[i] for i in range(self.n)]

    def dresid(self):
        Ay = [1]*self.n
        dmxtmulv(self.m, self.n, self.A, self.y, Ay)
        self.Hmul()
        self.rd = [self.cmod[i]*self.tau - self.s[i] - Ay[i]
                   for i in range(self.n)]

    def muresid(self):
        self.rmu = [self.mu - self.x[i]*self.s[i] for i in range(self.n)]
        if self.homo:
            self.hrmu = self.mu-self.kappa*self.tau
            self.rkxy = self.kappa+dot(self.cmod, self.x)-dot(self.b, self.y)

    def primal(self):
        self.Hmul()
        return sum([self.c[i]*self.x[i] + .5*self.Hx[i]*self.x[i] for i in range(self.n)])

    def dual(self):
        return sum([self.b[i]*self.y[i] for i in range(self.m)])-.5*sum([self.Hx[i]*self.x[i] for i in range(self.n)])

    def solvePC(self, gamma, corrector=0, dx0=[], ds0=[], dtau0=0, dkappa0=0):
        g1 = 1-gamma
        M = []
        if not self.tridiag:
            order = [0]*self.m
        else:
            Q = [0]*(self.m*self.m)
            DD = [0]*self.m
            EE = [0]*(self.m-1)
        if not corrector:
            self.Horder = [1]*self.n
            self.Hfact = [i for i in self.H]
            for i in range(self.n):
                self.Hfact[int(i*(i+3)/2)] += self.s[i]/self.x[i]
            if dsptrf('U', self.n, self.Hfact, self.Horder) != 0:
                raise 'Hfact is no good'
            for con in range(self.m):
                a = [self.A[con+self.m*k] for k in range(self.n)]
                dsptrs('U', self.n, 1, self.Hfact, self.Horder, a, self.n)
                for j in range(con+1):
                    M.append(dot(a, [self.A[j+self.m*k]
                             for k in range(self.n)]))
            print('M', M)
            back = 0
            if not self.tridiag:
                if self.m != 1:
                    back = dsptrf('U', self.m, M, order)  # factorise
            else:
                if self.m != 1:
                    back = trifact(self.m, M, Q, DD, EE)  # factorise
            if back != 0:
                raise 'Equations have gone singular %d' % back
            if not self.tridiag:
                self.M = M
                self.order = order
            else:
                self.Q = Q
                self.DD = DD
                self.EE = EE
        else:
            if not self.tridiag:
                M = self.M
                order = self.order
            else:
                Q = self.Q
                DD = self.DD
                EE = self.EE
        if not corrector:
            w1 = [self.rd[i]*g1 - (self.rmu[i]-g1*self.mu)/self.x[i]
                  for i in range(self.n)]
        else:
            w1 = [self.rd[i]*g1 - (self.rmu[i]-g1*self.mu -
                                   dx0[i]*ds0[i])/self.x[i] for i in range(self.n)]
        dsptrs('U', self.n, 1, self.Hfact, self.Horder, w1, self.n)
        dy = [0]*self.m
        dmxtmultv(self.n, self.m, self.A, w1, dy)
        dy = [dy[i]+self.rp[i]*g1 for i in range(self.m)]  # rhs
        if not self.tridiag:
            if self.m == 1:
                dy[0] /= M[0]
            else:
                dsptrs('U', self.m, 1, M, order, dy, self.m)  # solve for dy
        else:
            if self.m == 1:
                dy[0] /= M[0]
            else:
                trisolve(self.m, Q, DD, EE, dy)  # solve for dy
        if not self.homo:
            dx = [0]*self.n
            dmxtmulv(self.m, self.n, self.A, dy, dx)
            dsptrs('U', self.n, 1, self.Hfact, self.Horder, dx, self.n)
            dx = [dx[i]-w1[i] for i in range(self.n)]
            if not corrector:
                ds = [(self.rmu[i]-g1*self.mu-dx[i]*self.s[i])/self.x[i]
                      for i in range(self.n)]
            else:
                ds = [(self.rmu[i]-g1*self.mu-dx[i]*self.s[i]-dx0[i]
                       * ds0[i])/self.x[i] for i in range(self.n)]
            return(dx, ds, dy)
        else:
            db = [0]*self.m
            cx = [self.cmod[i] for i in range(self.n)]
            dsptrs('U', self.n, 1, self.Hfact, self.Horder, cx, self.n)
            dmxtmultv(self.n, self.m, self.A, cx, db)
            db = [db[i]+self.b[i] for i in range(self.m)]
            if not self.tridiag:
                if self.m == 1:
                    db[0] /= M[0]
                else:
                    dsptrs('U', self.m, 1, M, order, db, self.m)
            else:
                if self.m == 1:
                    db[0] /= M[0]
                else:
                    trisolve(self.m, Q, DD, EE, db)
            dx = [0]*self.n
            dmxtmulv(self.m, self.n, self.A, dy, dx)
            dsptrs('U', self.n, 1, self.Hfact, self.Horder, dx, self.n)
            dx = [dx[i]-w1[i] for i in range(self.n)]
            dc = [0]*self.n
            dmxtmulv(self.m, self.n, self.A, db, dc)
            dsptrs('U', self.n, 1, self.Hfact, self.Horder, dc, self.n)
            dc = [dc[i]-cx[i] for i in range(self.n)]
            cdx = dot(self.cmod, dx)
            bdy = dot(self.b, dy)
            cdc = dot(self.cmod, dc)
            bdb = dot(self.b, db)
            # -(cdx + dtaucdc) + bdy + dtaubdb = kappa + cx-by +(mu-taukappa - kappadtau)/tau
            if not corrector:
                dtau = (cdx-bdy+self.rkxy*g1+(self.hrmu-g1*self.mu) /
                        self.tau)/(bdb-cdc+self.kappa/self.tau)
            else:
                dtau = (cdx-bdy+self.rkxy*g1+(self.hrmu-g1*self.mu -
                        dtau0*dkappa0)/self.tau)/(bdb-cdc+self.kappa/self.tau)
            dx = [dx[i]+dtau*dc[i] for i in range(self.n)]
            dy = [dy[i]+dtau*db[i] for i in range(self.m)]
            if not corrector:
                ds = [(self.rmu[i]-g1*self.mu-dx[i]*self.s[i])/self.x[i]
                      for i in range(self.n)]
            else:
                ds = [(self.rmu[i]-g1*self.mu-dx[i]*self.s[i]-dx0[i]
                       * ds0[i])/self.x[i] for i in range(self.n)]
            if not corrector:
                dkappa = (self.hrmu-g1*self.mu-self.kappa*dtau)/self.tau
            else:
                dkappa = (self.hrmu-g1*self.mu-self.kappa *
                          dtau-dtau0*dkappa0)/self.tau
            return(dx, ds, dy, dtau, dkappa)


def IPopt_z(n, m, c, A, b, w, H=[], homogenous=1, sign=[], maxiter=100):
    """
    n variables, m linear constraints    
    w and c are vectors of length n
    A is m by n constraint matrix   A[j+m*i] ith element of jth constraint
    b is a vector of length m
    H is symmetric n by n           H is triangular, length is n*(n+1)/2

    Minimise c.x + 0.5*w.H.w
    subject to

    w[i]*sign[i]>=0     
    A.w=b


    To do linear programming just put H=[]
    Two methods;
    1) standard primal-dual interior point method using predictor-corrector
    2) homogenous primal-dual interior point method using predictor-corrector (default)

    sign[i] is either 1 or -1, to enable constraining w[i] positive or negative    
    """
    if H != []:
        Opt = QP(n, m, c, A, b, H)
    else:
        Opt = LP(n, m, c, A, b)
    if sign != []:
        Opt.sign = sign
        for i in range(n):
            if Opt.sign[i] == -1:
                Opt.x[i] = -1
                Opt.s[i] = -1
    Opt.homo = homogenous
    i = 0
    Opt.calcmu()
    mu0 = Opt.mu
    Opt.presid()
    Opt.dresid()
    Opt.muresid()
    Opt.results()
    rp0 = norm(Opt.rp)
    rd0 = norm(Opt.rd)
    compnow = comp0 = Opt.compl()
    # """
    toosmall = 0
    while 1:
        rp1 = norm(Opt.rp)/(1+rp0)
        rd1 = norm(Opt.rd)/(1+rd0)
        comp1 = Opt.compl()/(1+comp0)
        print('Iteration', i, rp1, rd1, comp1)
        if rp1 < sqeps and rd1 < sqeps and comp1 < sqeps:
            break
        if i > maxiter:
            break
        if Opt.homo:
            (dx1, ds1, dy1, dtau1, dkappa1) = Opt.solvePC(0)
        else:
            (dx1, ds1, dy1) = Opt.solvePC(0)
        if Opt.homo:
            alpha1 = .99*min(Opt.maxstep(dx1, ds1, dtau1, dkappa1))
        else:
            alpha1 = .99*min(Opt.maxstep(dx1, ds1))
        gamma = gfunc(alpha1)
        # gamma=1-min(comp1/compnow,1)
        if Opt.homo:
            (dx2, ds2, dy2, dtau2, dkappa2) = Opt.solvePC(gamma, 1, [alpha1*dx1[a] for a in range(n)], [alpha1*a for a in ds1],
                                                          alpha1*dtau1, alpha1*dkappa1)
        else:
            (dx2, ds2, dy2) = Opt.solvePC(gamma, 1, [
                alpha1*dx1[a] for a in range(n)], [alpha1*a for a in ds1])
        if Opt.homo:
            alpha2 = .99*min(Opt.maxstep(dx2, ds2, dtau2, dkappa2))
        else:
            alpha2 = .99*min(Opt.maxstep(dx2, ds2))
        if alpha1 > alpha2:
            print('Predictor')
            if not Opt.homo:
                Opt.update(dx1, dy1, ds1, alpha1)
            else:
                Opt.update(dx1, dy1, ds1, alpha1, dtau1, dkappa1)
        elif alpha2 > 2e-1:
            print('Corrector')
            if not Opt.homo:
                Opt.update(dx2, dy2, ds2, alpha2)
            else:
                Opt.update(dx2, dy2, ds2, alpha2, dtau2, dkappa2)
        else:
            print('Separate steps')
            if Opt.homo:
                (steph, steps, stepx) = Opt.maxstep(dx2, ds2, dtau2, dkappa2)
            else:
                (steph, steps, stepx) = Opt.maxstep(dx2, ds2)
            if Opt.homo:
                Opt.update_sep(dx2, dy2, ds2, .90*stepx, 1, .99 *
                               steps, dtau2, dkappa2, .99*steph)
            else:
                Opt.update_sep(dx2, dy2, ds2, .99*stepx,
                               1, .99*steps, 0, 0, .99*steph)
        Opt.calcmu()
        Opt.presid()
        Opt.dresid()
        Opt.muresid()
        Opt.results()
        compnow = comp1
        i += 1
    back = 1
    if rp1 < sqeps and rd1 < sqeps and comp1 < sqeps:
        print('Success')
        back = 0
        if Opt.homo:
            if Opt.tau < Opt.kappa and 1-Opt.tau/Opt.kappa > .9:
                print(1-Opt.tau/Opt.kappa)
                print('Infeasible')
                back = 6
            else:
                if not hasattr(Opt, 'cmod'):
                    Opt.cmod = Opt.c
                Opt.x = [Opt.x[i]/Opt.tau for i in range(Opt.n)]
                Opt.s = [Opt.s[i]/Opt.tau for i in range(Opt.n)]
                Opt.y = [Opt.y[i]/Opt.tau for i in range(Opt.m)]
                """
                print '%20s %20s %20s %20s'%('x','s','c','cmod')
                for i in range(Opt.n):
                    print '%20.8e %20.8e %20.8e %20.8e'%(Opt.x[i],Opt.s[i],Opt.c[i],Opt.cmod[i])
                print '%20s %20s'%('y','b')
                for i in range(Opt.m):
                    print '%20.8e %20.8e'%(Opt.y[i],Opt.b[i])
                """
        else:
            if not hasattr(Opt, 'cmod'):
                Opt.cmod = Opt.c
            """
            print '%20s %20s %20s %20s'%('x','s','c','cmod')
            for i in range(Opt.n):
                print '%20.8e %20.8e %20.8e %20.8e'%(Opt.x[i],Opt.s[i],Opt.c[i],Opt.cmod[i])
            print '%20s %20s'%('y','b')
            for i in range(Opt.m):
                print '%20.8e %20.8e'%(Opt.y[i],Opt.b[i])
            """
    else:
        print('Failed')
    if not back:
        for i in range(n):
            w[i] = Opt.x[i]
    return back


def IPopt(n, m, c, A, b, w, H=[], homogenous=1, sign=[], maxiter=100, L=[]):
    """Allow non-zero upper or lower bound (depends on sign list)"""
    if L == []:
        return IPopt_z(n, m, c, A, b, w, H, homogenous, sign, maxiter)
    CL = []
    if H != []:
        Sym_mult(n, H, L, CL)
        cc = [c[i]+CL[i] for i in range(n)]
    else:
        cc = c
    AL = [0]*m
    dmxtmultv(n, m, A, L, AL)
    bb = [b[i]-AL[i] for i in range(m)]
    back = IPopt_z(n, m, cc, A, bb, w, H, homogenous, sign, maxiter)
    for i in range(n):
        w[i] += L[i]
    return back


def SeqQP(n, m, c, A, b, w, H=[], homogenous=1, sign=[], maxiter=100, L=[]):
    """Try to do quadratic programming by sequential linear programming"""
    if H == []:
        return IPopt(n, m, c, A, b, w, H, homogenous, sign, maxiter, L)
    W = [0]*n
    if L == []:
        L = [0]*n
    if sign == []:
        sign = [1]*n
    cc = []
    Sym_mult(n, H, L, cc)
    cc = [c[i]+cc[i] for i in range(n)]
    back = IPopt(n, m, cc, A, b, W, [], homogenous, sign, maxiter, L)
    normw = 1
    ic = 0
    U0 = 1
    while ic <= 1 or (ic < 500 and back == 0):
        Sym_mult(n, H, W, cc)
        # print W
        U1 = dot(c, W)+0.5*dot(cc, W)
        print('Utility', U1)
        normw = differ(U0, U1)
        if normw < sqeps:
            break
        U0 = U1
        bb = [0]*m
        dmxtmultv(n, m, A, W, bb)
        w0 = [i for i in W]
        LL = [L[i]-W[i] for i in range(n)]
        bb = [b[i]-bb[i] for i in range(m)]
        cc = [c[i]+cc[i] for i in range(n)]
        back = IPopt(n, m, cc, A, bb, W, [], homogenous, sign, maxiter, LL)
        Sym_mult(n, H, W, bb)
        l = -dot(cc, W)/dot(W, bb)*.5  # Newton step along approximate path
        for i in range(n):
            if w0[i]+l*W[i] < L[i] and sign[i] == 1:
                raise '%20.8e %20.8e %20.8e %20.8e' % (
                    l, (L[i]-w0[i])/W[i], w0[i]+l*W[i], L[i])
                l = (L[i]-w0[i])/W[i]
            elif w0[i]+l*W[i] > L[i] and sign[i] == -1:
                raise '%20.8e %20.8e %20.8e %20.8e' % (
                    l, (L[i]-w0[i])/W[i], w0[i]+l*W[i], L[i])
                l = (L[i]-w0[i])/W[i]
        W = [w0[i]+l*W[i] for i in range(n)]
        print(ic, 'convergence', normw)
        ic += 1
    print('%d sequential iterations, convergence %20.8e' % (ic, normw))
    for i in range(n):
        w[i] = W[i]
    return back


if __name__ == '__main__':
    # Do sum(x)=-1 for x[i]<1e-1 i>0 x[0]<=-.5
    n = 1000
    c = [.01*(gfunc(i+1)+2) for i in range(n)]
    bench = [1.0/n]*n
    cextra = [0]*n
    m = 1
    A = [1]*n
    b = [1]
    H = [1e-2]*int(n*(n+1)/2)
    for i in range(n):
        H[int(i*(i+3)/2)] = .1*(1+i)
    gamma = 0
    Sym_mult(n, H, bench, cextra)
    for i in range(n):c[i] = -gamma/(1-gamma)*c[i]-cextra[i]
    w = [0]*n
    # [-0.5]+[3e-2]*(n-1)  # if sign[i] is 1 lower bound, -1 upper bound
    UL = []
#    SeqQP(n,m,c,A,b,w,H=H,homogenous=0,maxiter=100,sign=[-1]+[-1]*(n-1),L=UL)
    IPopt(n, m, c, A, b, w, H=H, homogenous=0,
          maxiter=100, sign=[1]+[1]*(n-1), L=UL)
    print(sum(w))
    if H != []:
        Hw = []
        Sym_mult(n, H, w, Hw)
        print(dot(c, w)+0.5*dot(Hw, w))
        for i in range(n):
            print('%20.8e %20.8e' % (w[i], c[i]+Hw[i]))
        Aw = [0]*m
        dmxtmultv(n, m, A, w, Aw)
        print(Aw, b)
