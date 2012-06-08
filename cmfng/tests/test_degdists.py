#!/usr/bin/env python
# coding: utf-8

"""Test for calculating degree distributions.
"""

import mfng
import unittest
import numpy
from known_values import known_values

class EnergyTests(unittest.TestCase):

    def testDegdistValues(self):
        "The calculated degree distribution should be correct."
        for i in range(3):
            divs   = known_values[i]["divs"]
            probs  = known_values[i]["probs"]
            result = known_values[i]["result"]

            pm = mfng.ProbMeasure(divs, probs)
            lpm = pm.iterate(4)
            ddnumpy = lpm.degdist_numpy(n=2000, maxdeg=20)
            dditerated = pm.degdist_iterated(n=2000, maxdeg=20)

            for i in range(20):
                self.assertAlmostEqual(ddnumpy[i]/result[i], 1, places=5)
            for i in range(20):
                self.assertAlmostEqual(dditerated[i]/result[i], 1, places=5)

def suite():
    degdist_suite = unittest.makeSuite(EnergyTests)
    return unittest.TestSuite([degdist_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())

if __name__ == "__main__":
    test()

