
import os
os.add_dll_directory('c:\\Users\\colin\\safeqp64')
from re import *
from BITAOPT import *
print(version())


class DATA:
    def __init__(self, ff):
        keyw = 0
        lastkey = 0
        self.Q = []
        for line in ff.readlines():
            if line.find('--------------') == 0:
                print('\x1b[1;1;32mBreak at %s \x1b[0;m' % line)
                break
            if not keyw:
                keyw = line.strip()
                # This allows numerical data across more than 1 line
                if search('^[0-9-]', keyw[0]):
                    print('\x1b[1;1;32mExtra line beginning with number for ',
                          '\x1b[1;1;31m', lastkey, '\x1b[1;1;36m', keyw.split()[0], '\x1b[0;m')
                    up = getattr(self, lastkey)
                    up += keyw.split()
                    setattr(self, lastkey, up)
                    keyw = 0
            else:
                setattr(self, keyw, line.strip().split())
                lastkey = keyw
                keyw = 0


ff = open('glog')  # Input all data from file glog
d = DATA(ff)
d.logfile = ''
for key in d.__dict__.keys():
    if key in 'tlen log n m nfac basket trades costs revise ls full round ncomp npiece nabs mabs longbasket shortbasket'.split():
        setattr(d, key, int(getattr(d, key)[0]))
    elif key in 'llambda gpower lpower C delta gamma Rmin Rmax LSValue LSValuel minRisk maxRisk'.split():
        setattr(d, key, float(getattr(d, key)[0]))
    elif key in 'min_holding min_trade DATA R FC FL SV mask L U A alpha bench Q initial min_lot size_lot Composites hpiece pgrad A_abs Abs_U Abs_L'.split():
        setattr(d, key, [float(i) for i in getattr(d, key)])
    elif key in 'I_A'.split():
        setattr(d, key, [int(i) for i in getattr(d, key)])

if d.nfac == -1:
    d.FC = []
    d.SV = []
    d.FL = []

w = []
gammaback = []
back = 0
if hasattr(d, 'C'):  # Gain Loss Variance
    back = GLOptimiseR(
        d.n,  # Number of assets
        d.nfac,  # Number of factors (-1 for full covariance)
        d.names,  # Asset Names
        # Number of dates for each asset (must be the same for each asset)
        d.tlen,
        # Historic return data matrix tlen by n, DATA[t+i*n] is t'th period for i'th asset
        d.DATA,
        d.R,  # Target return (Gain is Return-R, Loss is R-Return)
        d.C,  # Gain Loss utility is GAIN-C*LOSS
        d.llambda,  # multiplier of Gain-C*Loss in utility
        d.gpower,  # power of Gain in utility (experimental best set it to 1)
        d.lpower,  # power of Loss in utility (experimental best set it to 1)
        w,  # OUTPUT weights
        d.m,  # Number of linear constraints
        # Constraint matrix m by n, A[i+j*n] is i'th constraint value for j'th asset
        d.A,
        d.L,  # lower bounds for assets then lower bounds for constraints
        d.U,  # upper bounds for assets then upper bounds for constraints
        d.alpha,  # expected returns for assets
        d.bench,  # benchmark weights
        # the covariance matrix if nfac is -1 (array length n*(n+1)/2; packed symmetric matrix)
        d.Q,
        # the multiplier for return in the mean variance utility function is -gamma/(1-gamma)
        d.gamma,
        d.initial,  # initial portfolio weights
        d.delta,  # desired turnover (turnover is 0.5 sum|initial-w|
        d.basket,  # maximum desired number of non-zero asset weights
        d.trades,  # maximum desired number of non-zero asset trades
        d.revise,  # if revise is 1 we're doing a revision optimisation and initial must not be empty
        d.min_holding,  # smallest non-zero asset weight, empty means don't use
        d.min_trade,  # smallest non-zero trade, empty means don't use
        d.ls,  # if ls is 0 long only, if ls is 1 long short with long defining portfolio value, if ls is 2 long-short defines portfolio value
        d.full,  # if ls>0 upper portfolio value must be met
        d.Rmin,  # if ls>0 minimum short/long
        d.Rmax,  # if ls>0 maximum short/long
        d.LSValue,  # if ls>0 upper bound for portfolio value
        d.nabs,  # if ls>0 number of absolute constraints
        d.Abs_A,  # if ls>0 absolute constraint matrix for each asset
        d.mabs,  # number of linear constraints in A to be made absolute
        d.I_A,  # integer array of constraint numbers in A to be made absolute
        d.Abs_U,  # upper bounds for absolute constraints (length nabs+mabs)
        # if nfac!=-1 factor covariances array length nfac*(nfac+1)/2 symmtric packed
        d.FC,
        # matrix of factor loadings(betas) matrix n by nfac, FL[i+j*nfac] i'th beta for j'th asset
        d.FL,
        d.SV,  # array of specific variances
        d.mask,  # array length n of 0s and 1s which is applied to initial to enable assets to be excluded from turnover constraint
        d.log,  # control log output
        d.logfile,  # name of logfile
        d.longbasket,  # if ls >0 maximum number of positive non-zero weights
        d.shortbasket,  # if ls>0 maximum number of negative non-zero weights
        d.LSValuel,  # if ls>0 lower bound for portfolio value
        d.Abs_L,  # lower bounds for absolute constraints (length nabs+mabs)
        d.minRisk,  # if positive minimum desired risk
        d.maxRisk,  # if positive maximum desired risk
        gammaback  # OUTPUT the gamma used to get the correct risk
    )
else:  # Mean Variance Loss
    back = MVLOptimiseR(
        d.n,  # Number of assets
        d.nfac,  # Number of factors (-1 for full covariance)
        d.names,  # Asset Names
        # Number of dates for each asset (must be the same for each asset)
        d.tlen,
        # Historic return data matrix tlen by n, DATA[t+i*n] is t'th period for i'th asset
        d.DATA,
        d.R,  # Target return (Gain is Return-R, Loss is R-Return)
        d.llambda,  # multiplier of Gain-C*Loss in utility
        d.lpower,  # power of Loss in utility (experimental best set it to 1)
        w,  # OUTPUT weights
        d.m,  # Number of linear constraints
        # Constraint matrix m by n, A[i+j*n] is i'th constraint value for j'th asset
        d.A,
        d.L,  # lower bounds for assets then lower bounds for constraints
        d.U,  # upper bounds for assets then upper bounds for constraints
        d.alpha,  # expected returns for assets
        d.bench,  # benchmark weights
        # the covariance matrix if nfac is -1 (array length n*(n+1)/2; packed symmetric matrix)
        d.Q,
        # the multiplier for return in the mean variance utility function is -gamma/(1-gamma)
        d.gamma,
        d.initial,  # initial portfolio weights
        d.delta,  # desired turnover (turnover is 0.5 sum|initial-w|
        d.basket,  # maximum desired number of non-zero asset weights
        d.trades,  # maximum desired number of non-zero asset trades
        d.revise,  # if revise is 1 we're doing a revision optimisation and initial must not be empty
        d.min_holding,  # smallest non-zero asset weight, empty means don't use
        d.min_trade,  # smallest non-zero trade, empty means don't use
        d.ls,  # if ls is 0 long only, if ls is 1 long short with long defining portfolio value, if ls is 2 long-short defines portfolio value
        d.full,  # if ls>0 upper portfolio value must be met
        d.Rmin,  # if ls>0 minimum short/long
        d.Rmax,  # if ls>0 maximum short/long
        d.LSValue,  # if ls>0 upper bound for portfolio value
        d.nabs,  # if ls>0 number of absolute constraints
        d.Abs_A,  # if ls>0 absolute constraint matrix for each asset
        d.mabs,  # number of linear constraints in A to be made absolute
        d.I_A,  # integer array of constraint numbers in A to be made absolute
        d.Abs_U,  # upper bounds for absolute constraints (length nabs+mabs)
        # if nfac!=-1 factor covariances array length nfac*(nfac+1)/2 symmtric packed
        d.FC,
        # matrix of factor loadings(betas) matrix n by nfac, FL[i+j*nfac] i'th beta for j'th asset
        d.FL,
        d.SV,  # array of specific variances
        d.mask,  # array length n of 0s and 1s which is applied to initial to enable assets to be excluded from turnover constraint
        d.log,  # control log output
        d.logfile,  # name of logfile
        d.longbasket,  # if ls >0 maximum number of positive non-zero weights
        d.shortbasket,  # if ls>0 maximum number of negative non-zero weights
        d.LSValuel,  # if ls>0 lower bound for portfolio value
        d.Abs_L,  # lower bounds for absolute constraints (length nabs+mabs)
        d.minRisk,  # if positive minimum desired risk
        d.maxRisk,  # if positive maximum desired risk
        gammaback  # OUTPUT the gamma used to get the correct risk
    )

print(Return_Message(back))
prob = []
gain = []
loss = []
marggain = []
margloss = []

GLProp(d.n, d.names, d.tlen, d.DATA, d.R, d.gpower,
       d.lpower, prob, gain, loss, w, marggain, margloss)

print('Probability of gain', prob[0])
print('Expected gain', gain[0])
print('Expected loss', loss[0])

print(('%20s %20s %20s %20s' % ('Asset', 'weight', 'marg gain', 'marg loss')))
for i in range(d.n):
    print(('%20s %20.5f %20.5f %20.5f' %
          (d.names[i], w[i], marggain[i], margloss[i])))
