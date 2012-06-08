#!/usr/bin/env python
# coding: utf-8

"""Unittest for analyzer.py."""

import cmfng.analyzer as an
import numpy
import unittest
import os

class LablesTests(unittest.TestCase):

    project = "test"

    def testLoad(self):
        r = an.Runs(self.project)
        labels = r.set_labels()
        self.assertEqual(len(labels), 2)
        self.assertEqual(labels[0], "2000_001")
        self.assertEqual(labels[1], "2000_002")
        r.loglog1()
        r.loglog()
        r.properties(10)
        r.degdist()

    def testNoLabelsError(self):
        r = an.Runs(self.project)
        self.assertRaises(an.NoLabelsError, r.loglog1)
        self.assertRaises(an.NoLabelsError, r.loglog)
        self.assertRaises(an.NoLabelsError, r.properties)
        self.assertRaises(an.NoLabelsError, r.degdist)
        r.loglog1(label="2000_001")
        r.degdist(label="2000_001")

if __name__ == "__main__":
    actual_dir = os.path.abspath(".")
    test_dir = os.path.split(__file__)[0]
    if test_dir:
        os.chdir(test_dir)
    unittest.main()
    os.chdir(actual_dir)

