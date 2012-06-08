#!/usr/bin/env python
# coding: utf-8

"""Unittest for the Generator class."""

from __future__ import division
from __future__ import print_function
import cmfng
import unittest
import numpy
import random
import shutil
import os
from cmfng import analyzer
from known_values import known_values

class GeneratorTests(unittest.TestCase):
    def testGenerationandRun(self):
        "The Generator should work with all the properties and create labels/files."
        dir_ = "project_base"
        if os.path.isdir(dir_):
            shutil.rmtree(dir_)
        properties = (
            cmfng.AverageDegree(value=10),
            cmfng.DistributionFunction("k**-2", maxdeg=10, mindeg=1),
            cmfng.LogBinnedDegDist(degrees = [1] * 5 + [12]),
            )
        for property_ in properties:
            generator = cmfng.Generator(T0=.2, steps=10, Tlimit=.04)
            generator.append_property(property_)
            generator.go()
        for m in [2, 3, 4]:
            generator = cmfng.Generator(T0=.2, steps=10, Tlimit=.04, m=m)
            generator.append_property(property_)
            generator.go()
        generator = cmfng.Generator(T0=.2, steps=10, Tlimit=.04, n=400)
        generator.append_property(properties[0])
        generator.go()
        runs = analyzer.Runs()
        self.assertTrue(os.path.isfile("{0}/runs.py".format(dir_)))
        labels = runs.set_labels()
        self.assertEqual(len(labels), 7)
        self.assertEqual(labels, ["200_{0:03}".format(i) for i in range(1,7)]+["400_001"])

    def testAccept(self):
        "Generator should accept better ProbMeasures."
        generator = cmfng.Generator(T0=.2, steps=10, Tlimit=.04, m=3)
        generator.append_property(cmfng.DistributionFunction("k**-2", maxdeg=10, mindeg=1))
        for i in range(3):
            divs   = known_values[i]["divs"]
            probs  = known_values[i]["probs"]
            energy = known_values[i]["energy"]
            self.assertTrue(generator.has_accepted(cmfng.ProbMeasure(divs, probs)))

    def testBadM(self):
        "m should be integer and greater then 1"
        bad_m_values = [2.0, 2.3, 1, 0, -3, ]
        for m in bad_m_values:
            self.assertRaises(AssertionError, cmfng.Generator,
                    T0=.2, steps=10, Tlimit=.04, m=m)

    def testBadDivexponents(self):
        "divexponent should be positive integer"
        bad_divexponents = [1.0, 1.2, -1, 0]
        for divexponent in bad_divexponents:
            self.assertRaises(AssertionError, cmfng.Generator, T0=.2, steps=10, Tlimit=.04, divexponent=divexponent)

def suite():
    generaror_suite = unittest.makeSuite(GeneratorTests)
    return unittest.TestSuite([generaror_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())

if __name__ == "__main__":
    test()

