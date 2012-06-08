import unittest
from cmfng.tests import (test_analyzer, test_degdists,
        test_DistributionFuntion, test_generation, test_logbinned_degdist,
        test_probmeasure, test_step_number,)

def suite():
    return unittest.TestSuite([
        test_degdists.suite(),
        test_generation.suite(),
        test_DistributionFuntion.suite(),
        test_logbinned_degdist.suite(),
        test_probmeasure.suite(),
        test_step_number.suite(),
        ])

def run_tests(verbosity=1):
    try:
        # Support for testoob to have nice colored output
        import testoob
        testoob.main(suite())
    except ImportError:
        runner = unittest.TextTestRunner(verbosity=verbosity)
        runner.run(suite())

# Make nosetest skip run_tests
run_tests.__test__ = False

if __name__ == "__main__":
    run_tests(255)
    print("Run test_analyzer.py separately.")

