#!/usr/bin/env python
# -*- coding: utf-8 -*-

try:
    from numpy import array
except ImportError:
    array = lambda x: x

runs={}
runs["2000_001"] = dict(
  # Settings
     T0=0.2,  Tfactor=0.7, Tlimit=2e-05,
     m=3, K=4,
     n=2000,
     divexponent=7,
     properties= [
        "DistributionFunction('k**-2', maxdeg=1999, mindeg=1, bigfloat=False, division=True)"
     ],
  # Results
     divs= [0.033743820791589288, 0.64256466270430945, 1.0],
     probs=\
       array([[ 0.11325525,  0.10898781,  0.08560732],
       [ 0.10898781,  0.11897457,  0.12681913],
       [ 0.08560732,  0.12681913,  0.12494168]]),
     #Energies
       Ei =1998.303391, #Initial
       Ef =1997.888507, #Final
     runtime =  0.550568, # minutes (2012.01.26 13:51:04 -- 2012.01.26 13:51:37)
     steps =  52,

     accept_reject = "AA AA Aa AA AA Aa Aa a. aa aA .a .a .a Aa .a A. AA Aa .. .. A. .A Aa .. .. Aa ",

)

runs["2000_002"] = dict(
  # Settings
     T0=0.2,  Tfactor=0.7, Tlimit=2e-05,
     m=3, K=4,
     n=2000,
     divexponent=7,
     properties= [
        "DistributionFunction('k**-2', maxdeg=1999, mindeg=1, bigfloat=False, division=True)"
     ],
  # Results
     divs= [0.23255815009953937, 0.84800719775645905, 1.0],
     probs=\
       array([[ 0.10086081,  0.12910521,  0.09643349],
       [ 0.12910521,  0.11846284,  0.11157292],
       [ 0.09643349,  0.11157292,  0.1064531 ]]),
     #Energies
       Ei =1998.303391, #Initial
       Ef =1998.082190, #Final
     runtime =  0.551340, # minutes (2012.01.26 13:51:37 -- 2012.01.26 13:52:10)
     steps =  52,

     accept_reject = "Aa Aa .a aA aA aa Aa .a .a A. .A AA A. .. Aa .A Aa Aa A. A. A. .a A. .. .a Aa ",

)

