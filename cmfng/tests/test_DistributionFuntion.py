#!/usr/bin/env python
# coding: utf-8

"""Unittest for the target property classes of mfng.py."""

import cmfng
import unittest
import numpy
from known_values import known_values

class CreationTests(unittest.TestCase):

    def testDistributionFunctionValues(self):
        "The value of the DistributionFunction should be correct."
        distfunc = cmfng.DistributionFunction("2*numpy.exp(-k)", 100)
        self.assertAlmostEqual(distfunc.value(numpy.log(8)), 0.25)

#    def testEnergy(self):
#        df=cmfng.DistributionFunction("2*numpy.exp(-k)", 100)
#        pm=cmfng.ProbMeasure()
#        lpm=pm.iterate(3)
#        self.assertAlmostEqual(df.energy(lpm, 400)[0], -97.436661566437408)
#
#    def testEnergy2(self):
#        divs= [0.13193312663039208, 1.0]
#        probs=numpy.array([[ 0.25904704,  0.30615741],
#               [ 0.30615741,  0.12863815]])
#        pm=cmfng.ProbMeasure(divs, probs)
#        lpm=pm.iterate(3)
#        df=cmfng.DistributionFunction("k**-3", 100, mindeg=1)
#        self.assertAlmostEqual(df.energy(lpm, 400)[0], -92.244977910442273)

class EnergyTests(unittest.TestCase):

    def testDistributionFunctionValues(self):
        "The energy of the DistributionFunction should be correct."
        for i in range(2):
            divs   = known_values[i]["divs"]
            probs  = known_values[i]["probs"]
            energy = known_values[i]["energy"]

            pm = cmfng.ProbMeasure(divs, probs)
            lpm = pm.iterate(4)
            prop = cmfng.DistributionFunction('k**-2', maxdeg=2000-1, mindeg=1)
            self.assertAlmostEqual(prop.energy(lpm, 2000)[0], energy)
            prop = cmfng.DistributionFunctionC('k**-2', maxdeg=2000-1, mindeg=1)
            self.assertAlmostEqual(prop.energy(pm, 2000)[0], energy, places=4)

def suite():
    energy_suite = unittest.makeSuite(EnergyTests)
    creation_suite = unittest.makeSuite(CreationTests)
    return unittest.TestSuite([energy_suite, creation_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())

if __name__ == "__main__":
    test()

