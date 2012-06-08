#!/usr/bin/env python
# coding: utf-8

"""Unittest for mfng.py."""

import cmfng
import numpy
import unittest

class Given_m_Tests(unittest.TestCase):
    m_max = 8

    def testMOne(self):
        self.assertRaises(AssertionError, cmfng.simpleProbMeasure, m=1)

    def testProbsShape(self):
        for m in range(2,self.m_max):
            pm = cmfng.simpleProbMeasure(m=m)
            shape = pm.probs.shape
            self.assertEqual(shape, (m,m))

    def testProbsValues(self):
        for m in range(2,self.m_max):
            prob = 1./m**2
            for i in range(m):
                for j in range(m):
                    pm = cmfng.simpleProbMeasure(m=m)
                    p = pm.probs[i,j]
                    self.assertAlmostEqual(p, prob)

    def testDivsLength(self):
        for m in range(2,self.m_max):
            pm = cmfng.simpleProbMeasure(m=m)
            length = len(pm.divs)
            self.assertEqual(length, m)

    def testDivsValues(self):
        for m in range(2,self.m_max):
            pm = cmfng.simpleProbMeasure(m=m)
            divs = pm.divs
            value = 0.
            delta = 1./m
            for i in range(m):
                value += delta
                self.assertAlmostEqual(value, divs[i])


class IterationValuesTests(unittest.TestCase):
    m_max = 8
    K_max = 4
    known_values = (
        dict(
            m=2,
            K=2,
            divs =  [0.25, 0.5, 0.75, 1],
            probs = numpy.array(
                [
                 [1./16, 1./16, 1./16, 1./16, 1./16, ],
                 [1./16, 1./16, 1./16, 1./16, 1./16, ],
                 [1./16, 1./16, 1./16, 1./16, 1./16, ],
                 [1./16, 1./16, 1./16, 1./16, 1./16, ],
                ]
                )
        ),
        )

    def testTooManyDivs(self):
        m  = 4
        pm = cmfng.simpleProbMeasure(m=m)
        K = 4
        # 4**4 = 216
        self.assertRaises(ValueError, pm.iterate, K=K)

    def testMaxDivs(self):
        for m in range(2,self.m_max):
            pm = cmfng.simpleProbMeasure(m=m)
            for K in range(1,self.K_max):
                self.assertRaises(ValueError, pm.iterate, K=K, maxdivs=m**K-1)

    def testProbsShape(self):
        for m in range(2,self.m_max):
            pm = cmfng.simpleProbMeasure(m=m)
            for K in range(1,self.K_max):
                ipm = pm.iterate(K=K, maxdivs=4100)
                shape = ipm.probs.shape
                self.assertEqual(shape, (m**K,m**K))

    def testDivsLength(self):
        for m in range(2,self.m_max):
            pm = cmfng.simpleProbMeasure(m=m)
            for K in range(1,self.K_max):
                ipm = pm.iterate(K=K, maxdivs=4100)
                length = len(ipm.divs)
                self.assertEqual(length, m**K)

    def testIteratedValue(self):
        for value in self.known_values:
            m = value["m"]
            K = value["K"]
            divs = value["divs"]
            pm = cmfng.simpleProbMeasure(m=m)
            ipm = pm.iterate(K=K)
            for i, j in zip(ipm.divs, divs):
                self.assertAlmostEqual(i, j)

class AvgDegreeTests(unittest.TestCase):
    def testAvgDegree(self):
        pm = cmfng.simpleProbMeasure(m=3)
        self.assertAlmostEqual(pm.avg_degree(9001), 1000)

class ProbsDivsTests(unittest.TestCase):
    #TODO "An older test."
    def setUp(self):
        self.pm = cmfng.ProbMeasure(
            numpy.array([.3, 1]),
            numpy.array([[0.1 , .25],
                         [0.25, 0.4]])
            )
        self.iterated_pm = self.pm.iterate(2)

    def testProbs(self):
        diff = numpy.sum(self.pm.probs-
        numpy.array([[ 0.1  ,  0.25],
               [ 0.25 ,  0.4 ]]))
        self.assertTrue(diff < 1e7)

    def testDivs(self):
        diff = numpy.sum(numpy.array(self.pm.divs)-
          numpy.array([0.3, 1]))
        self.assertTrue(diff < 1e7)

    def testIteratedProbs(self):
        diff = numpy.sum(self.iterated_pm.probs-
        numpy.array([[ 0.01  ,  0.025 ,  0.025 ,  0.0625],
                     [ 0.025 ,  0.04  ,  0.0625,  0.1   ],
                     [ 0.025 ,  0.0625,  0.04  ,  0.1   ],
                     [ 0.0625,  0.1   ,  0.1   ,  0.16  ]]))
        self.assertTrue(diff < 1e7)

    def testIteratedDivs(self):
        diff = numpy.sum(numpy.array(self.iterated_pm.divs)-
          numpy.array([0.09, 0.3, 0.51, 1]))
        self.assertTrue(diff < 1e7)


def suite():
    givenm_suite = unittest.makeSuite(Given_m_Tests)
    iteration_suite = unittest.makeSuite(IterationValuesTests)
    avgdeg_suite = unittest.makeSuite(AvgDegreeTests)
    probsdivs_suite = unittest.makeSuite(ProbsDivsTests)
    return unittest.TestSuite([givenm_suite, iteration_suite, avgdeg_suite,
        probsdivs_suite])

def test():
    runner = unittest.TextTestRunner()
    runner.run(suite())

if __name__ == "__main__":
    test()

