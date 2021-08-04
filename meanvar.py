
OPTPATH='c:\\Users\\colin\\safeqp64'
from os import add_dll_directory
add_dll_directory(OPTPATH)
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


ff = open('matlog')
d = DATA(ff)
d.logfile = ''
for key in d.__dict__.keys():
    if key in 'log n m nfac basket tradenum costs revise ls full round ncomp npiece nabs mabs longbasket shortbasket tradebuy tradesell'.split():
        setattr(d, key, int(getattr(d, key)[0]))
    elif key in 'delta kappa gamma min_holding min_trade rmin rmax value valuel minRisk maxRisk ShortCostScale'.split():
        setattr(d, key, float(getattr(d, key)[0]))
    elif key in 'DATA R C FC FL SV mask L U A alpha bench Q initial buy sell min_lot size_lot Composites hpiece pgrad A_abs Abs_U'.split():
        setattr(d, key, [float(i) for i in getattr(d, key)])
    elif key in 'I_A'.split():
        setattr(d, key, [int(i) for i in getattr(d, key)])
d.zetaS = 1
d.zetaF = 1
d.downrisk = 0
d.downfactor = 3
w = []
shake = []
ogamma = []

# Minimise utility which contains risk, return and cost subject to a basket constraint
# The outcome weights will go to list w
back = Optimise_internalCVPAFbl(
    d.n,  # Number of assets
    d.nfac,  # Number of factors (-1 for full covariance)
    d.names,  # Asset Names
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
    # the multiplier for return in the utility function is -gamma/(1-gamma)
    d.gamma,
    d.initial,  # initial portfolio weights
    d.delta,  # desired turnover (turnover is 0.5 sum|initial-w|
    d.buy,  # transaction buy costs for each asset (+ve)
    d.sell,  # transaction sell costs for each asset (+ve)
    # the multiplier for cost in the utility function is kappa/(1-kappa) if kappa is -1 set kappa to -gamma
    d.kappa,
    d.basket,  # maximum desired number of non-zero asset weights
    d.tradenum,  # maximum desried number of non-zero trades
    d.revise,  # if revise is 1 we're doing a revision optimisation and initial must not be empty
    d.costs,  # if costs is 0 ignore costs, if costs is 1 include costs, if costs is 2 include costs in  the budget constraint
    d.min_holding,  # smallest non-zero asset weight, empty means don't use
    d.min_trade,  # smallest non-zero trade, empty means don't use
    d.ls,  # if ls is 0 long only, if ls is 1 long short with long defining portfolio value, if ls is 2 long-short defines portfolio value
    d.full,  # if ls>0 upper portfolio value must be met
    d.rmin,  # if ls>0 minimum short/long
    d.rmax,  # if ls>0 maximum short/long
    d.round,  # use round lots if round is 1
    d.min_lot,  # lowest lot size for each asset
    d.size_lot,  # subsequent lot size for each asset
    shake,  # array, shake[i]=-1 if asset i was rounded correctly
    d.ncomp,  # number for assets which are composite
    d.Composites,  # matrix defining composite assets weights for each asset
    d.value,  # if ls>0 upper bound for portfolio value
    d.npiece,  # if npiece is 1 use piecewise costs
    d.hpiece,  # matrix for x-coordinates for piecewise costs for each asset npiece by n
    d.pgrad,  # matrix piecewise cost gradients for each asset npiece by n
    d.nabs,  # if ls>0 number of absolute constraints
    d.A_abs,  # if ls>0 absolute constraint matrix for each asset
    d.mabs,  # number of linear constraints in A to be made absolute
    d.I_A,  # integer array of constraint numbers in A to be made absolute
    d.Abs_U,  # upper bounds for absolute constraints (length nabs+mabs)
    # if nfac!=-1 factor covariances array length nfac*(nfac+1)/2 symmtric packed
    d.FC,
    # matrix of factor loadings(betas) matrix n by nfac, FL[i+j*nfac] i'th beta for j'th asset
    d.FL,
    d.SV,  # array of specific variances
    d.minRisk,  # if positive minimum desired risk
    d.maxRisk,  # if positive maximum desried risk, gamma is varied to get the desrired risk and optimal gamma is returned in ogamma
    ogamma,  # OUTPUT the gamma used to get the correct risk
    d.mask,  # array length n of 0s and 1s which is applied to initial to enable assets to be excluded from turnover constraint
    d.log,  # control log output
    d.logfile,  # name of logfile
    d.downrisk,  # experimental keep at 0
    d.downfactor,  # experimental keep at 3
    d.longbasket,  # if ls >0 maximum number of positive non-zero weights
    d.shortbasket,  # if ls>0 maximum number of negative non-zero weights
    d.tradebuy,  # maximum number of buys
    d.tradesell,  # maximum number of sells
    d.zetaS,  # experimental set to 1
    d.zetaF,  # experimental set to 1
    d.ShortCostScale,  # scale short costs by this (1 means same as long costs)
    d.valuel,  # if ls>0 lower bound for portfolio value
    d.Abs_L  # lower bounds for absolute constraints (length nabs+mabs)
)


print(Return_Message(back))
arisk = []
risk = []
Rrisk = []
brisk = []
pbeta = []

# Calculate some risks using the optimal weights
Get_RisksC(d.n, d.nfac, d.Q, w, d.bench, arisk, risk,
           Rrisk, brisk, pbeta, d.ncomp, d.Composites)

print('Portfolio Absolute Risk', arisk[0])
print('Portfolio Relative Risk', risk[0])
print('Portfolio Residual Risk', Rrisk[0])
print('Portfolio Benchmark Risk', brisk[0])
print('Portfolio beta', pbeta[0])
