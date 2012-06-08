#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
import mfng

for avg_degree in [2, 100]:
    T0 = 0.2
    Tlimit = T0/10000
    generator = mfng.Generator(T0=T0, steps=100, Tlimit=Tlimit,
            m=3, K=3,
            n=8000,
            divexponent = 7,
            project = "base",
            )
    generator.append_property(
            mfng.AverageDegree(avg_degree)
            )
    generator.go()

