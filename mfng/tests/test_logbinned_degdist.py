#!/usr/bin/env python
# coding: utf-8

"""Unittest for the LogBinnedDegDist property class of mfng.py."""

import cmfng
import unittest
import numpy
import random

class LogBinnedDegDistTests(unittest.TestCase):
    known_values = (
        ([1,1,2,4], [.5,.25,.25]),
        ([1,1,1,1,2,2,4,8], [.5,.25,.125,.125]),
        ([1,2,4,8,16], [1./5]*5),
        (range(0,16), [.0625,.125,.25,.5]),
        ([1,16], [.5,0,0,0,.5]),

        # With the same binned distribution.
        ([ 8,256], [0,0,0,.5,0,0,0,0,.5]),
        ([10,299], [0,0,0,.5,0,0,0,0,.5]),
        ([15,511], [0,0,0,.5,0,0,0,0,.5]),
        )
    def testBinnedDegdistValues(self):
        "LogBinnedDegDist should have the proper binned_degdist"
        for degrees, binned_degdist in self.known_values:
            random.shuffle(degrees)
            lbdd = cmfng.LogBinnedDegDist(degrees).binned_degdist
            self.assertEqual(len(lbdd), len(binned_degdist))
            for i in range(len(binned_degdist)):
                self.assertAlmostEqual(lbdd[i], binned_degdist[i])

class LogBinnerTests(unittest.TestCase):

    known_values = (
        ([1] * 5, [1, 2, 2]),
        ([1] * 7, [1, 2, 4]),
        ([1] * 8, [1, 2, 4, 1]),
        (range(8), [0, 3, 18, 7]),
        ([.5, .25, .125, .125], [.5, .375, .125]),
        ([0, 0, .5, .25, 0, .25], [0, .5, .5]),
        ([.2] * 5, [.2, .4, .4]),
        )

    def testLogBinnerValues(self):
        "log_binner should give the proper values."
        for list_, binned in self.known_values:
            self.assertEqual(cmfng.log_binner(list_), binned)

    bad_args = (
        [0,1,2,"a"],
        [1j],
        "alma",
        1,
        )
    def testLogBinnerBadValues(self):
        "log_binner should raise ValueError for these arguments."
        for value in self.bad_args:
            self.assertRaises(ValueError, cmfng.log_binner, value)

def suite():
    lbdegdist_suite = unittest.makeSuite(LogBinnedDegDistTests)
    logbinner_suite = unittest.makeSuite(LogBinnerTests)
    return unittest.TestSuite([logbinner_suite, lbdegdist_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())

if __name__ == "__main__":
    test()

