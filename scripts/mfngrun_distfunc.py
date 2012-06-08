#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
import cmfng
import sys

n = 2000
maxdeg = n - 1

for steps in [1000]:
    T0 = 0.2
    Tlimit = T0/10000
    generator = cmfng.Generator(T0=T0, steps=steps, Tlimit=Tlimit,
            m=3, K=4,
            n=n,
            divexponent=7,
            project="distfunc_division_experimental",
            )
    generator.append_property(
            cmfng.DistributionFunctionC(
                "k**-2",
                maxdeg=maxdeg, mindeg=1, K=generator.K
                )
            )
    generator.go()

