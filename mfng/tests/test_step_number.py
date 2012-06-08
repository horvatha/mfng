#!/usr/bin/env python
# coding: utf-8

"""Unittest for the Generator class."""

from __future__ import division
from __future__ import print_function
import cmfng
import unittest

class GenerationStepsTests(unittest.TestCase):

    known_values = ( # T0/Tlimit, Tfactor, steps
            (100, 0.987, 704),
            (100, 0.988, 764),
            (100, 0.989, 834),
            (100, 0.99, 918),
            (100, 0.991, 1020),
            (100, 0.992, 1148),
            (100, 0.993, 1312),
            (100, 0.994, 1532),
            (100, 0.995, 1838),
            (100, 0.996, 2298),
            (100, 0.997, 3066),
            (100, 0.998, 4602),
            (100, 0.999, 9206),
            (1000, 0.987, 1056),
            (1000, 0.988, 1146),
            (1000, 0.989, 1250),
            (1000, 0.99, 1376),
            (1000, 0.991, 1530),
            (1000, 0.992, 1722),
            (1000, 0.993, 1968),
            (1000, 0.994, 2296),
            (1000, 0.995, 2758),
            (1000, 0.996, 3448),
            (1000, 0.997, 4600),
            (1000, 0.998, 6902),
            (1000, 0.999, 13810),
            (10000, 0.987, 1408),
            (10000, 0.988, 1526),
            (10000, 0.989, 1666),
            (10000, 0.99, 1834),
            (10000, 0.991, 2038),
            (10000, 0.992, 2294),
            (10000, 0.993, 2624),
            (10000, 0.994, 3062),
            (10000, 0.995, 3676),
            (10000, 0.996, 4596),
            (10000, 0.997, 6132),
            (10000, 0.998, 9202),
            (10000, 0.999, 18412),
            (100000, 0.987, 1760),
            (100000, 0.988, 1908),
            (100000, 0.989, 2082),
            (100000, 0.99, 2292),
            (100000, 0.991, 2548),
            (100000, 0.992, 2868),
            (100000, 0.993, 3278),
            (100000, 0.994, 3828),
            (100000, 0.995, 4594),
            (100000, 0.996, 5746),
            (100000, 0.997, 7664),
            (100000, 0.998, 11502),
            (100000, 0.999, 23016),
            )

    def testStepNumber(self):
        "The cmfng.steps sould give proper values."
        for range_, Tfactor, steps in self.known_values:
            T0 = .2
            self.assertEqual(cmfng.steps(T0, Tfactor, T0/range_), steps)
            self.assertEqual(cmfng.steps(T0, Tfactor, T0/float(range_)), steps)

    def testTfactor(self):
        "The cmfng.Tfactor sould give proper values."
        for quotient in (1/10**_ for _ in range(3, 9)):
            for steps in (10**_ for _ in range(3, 9)):
                T0 = .2
                steps_ = cmfng.steps(T0, cmfng.get_Tfactor(quotient, steps), T0*quotient)
                if (steps-steps_) != 0:
                    print(quotient, steps_)
                self.assertEqual(steps_, steps)

def suite():
    generationsteps_suite = unittest.makeSuite(GenerationStepsTests)
    return unittest.TestSuite([generationsteps_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())

if __name__ == "__main__":
    test()


