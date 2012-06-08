#!/usr/bin/env python
# coding: utf-8
"""
Multifractal network generator.

See:
Palla - Vicsek - Lovász: Multifractal network generator, PNAS, 2010
"""

from __future__ import with_statement # Not in Python 2.5
from __future__ import print_function
from __future__ import division

from random import randrange, random
import os
import numpy
from bisect import bisect_left
try:
    import igraph
except ImportError:
    pass
import sys
import shelve
import json
import math
import re
import string

if False: # Analyzing memory usage
    from guppy import hpy
    heapy = hpy()
    from gc import collect
    global f
    f = open("memstat_out.txt", "a")
    import os
    memstat = "memstat -p %d" % os.getpid()
    def printmem(verbose=False, fi = f, text=""):
        stdout = os.popen(memstat)
        lines = stdout.readlines()
        if verbose:
            out = "".join(lines)
        else:
            lastline = lines[-1].split()
            out = "[%s%s(%s]" % (lines[0].split()[0], lastline[0], lastline[2])
        print(out)
        if fi:
            fi.write(text, out)

class Output:
    def __init__(self, files=None, with_stdout=True):
        if files is None:
            files = []
        elif isinstance(files, str):
            files = [open(files, "a")]
        elif isinstance(files, list):
            assert all(isinstance(f, file) for f in files)
        else:
            raise ValueError("files must be None, string or list")
        if with_stdout:
            files.append(sys.stdout)
        self.files = files

    def __del__(self):
        for f in self.files:
            if f is not sys.stdout:
                f.close()

    def write(self, text):
        """Write the text into the file objects."""
        for f in self.files:
            print(text, file=f)


class NoOutput:
    def write(self, text):
        pass

from time import localtime, time, strftime
hms = lambda _time: ".".join("%02d" % i for i in localtime(_time)[:3]) + \
      " " +  ":".join("%2.2d" % i for i in localtime(_time)[3:6])


class ProbMeasure(object):
    """Probability measure

    See:
    Palla - Vicsek - Lovász: Multifractal network generator, PNAS, 2010

    """
    def __init__(self, divs, probs, output=None, verbose=True, iterated=False):
        '''Creates the divs and matrix.

        divs: list or None, default None
           the divisors of the [0,1] interval. E.g.:
           [0.2, 0.5]
           or
           [0.2, 0.5, 1]

        probs: numpy.ndarray or string
           If string, the values of the matrix in the form of string like
           """
           .3
           .4,.6
           .2,.3,.1
           """
           __init__ creates and normalizes the symmetric matrix.

        m: positive integer
           If no divs and probs given, it generates a matrix with
           equal probabilities and divs with equal intervals.

        '''
        self.m = len(divs)
        self.divs = list(divs)
        # TODO assert isinstance(probs, numpy.ndarray)
        self.probs = probs
        self.verbose = verbose
        self.iterated = iterated
        self.probs /= numpy.sum(self.probs)
        if output is None:
            if verbose:
                self.out = Output()
            else:
                self.out = NoOutput()
        else:
            self.out = output

    def __str__(self):
        pm_type = "Link probability measure" if self.iterated else "Generator measure"
        return "{0} with {1}x{1} elements.".format(pm_type, self.m)

    def summary(self):
        pm_type = "Link probability measure" if self.iterated else "Generator measure"
        return "{0} with {1}x{1} elements.\ndivs={2}\nprobs={3}".format(
                    pm_type, self.m, self.divs, self.probs
                    )

    def avg_degree(self, n):
        """Returns the average degree.

        It should be use with the (iterated) link probability measure.

        Parameters:
          n: the size of the network (the number of vertices)

        Returns:
          The average degree.
        """

        assert isinstance(n, int)
        divs = self.divs
        n_intervals = len(divs)
        lengths = [divs[i] - divs[i-1] for i in xrange(1, n_intervals)]
        lengths.insert(0, divs[0])
        # Eq. 5, where ...
        avgdeg_i = [sum([self.probs[i][j] * lengths[j]
                  for j in xrange(n_intervals)]) for i in xrange(n_intervals)]
        avgdeg = (n-1)*sum(avgdeg_i[i] * lengths[i] for i in  xrange(n_intervals))
        return avgdeg

    def degdist_iterated(self, maxdeg, n, mindeg=0, K=4):
        """Returns with the degree distribution of the iterated prob. measure.

        It uses a fast external program.
        """

        divs_a = " ".join(map(str,self.divs))
        divs_a = "0.0 " + divs_a
        tmp_probs_a = []
        for i in range(0,len(self.probs)):
            tmp_probs_a.append( " ".join(map(str,self.probs[i])) )
        probs_a = "\n".join(map(str,tmp_probs_a))

        stdin_str = divs_a + "\n" + probs_a + "\n"

        fd=open("x.tmp","w")
        fd.write(stdin_str)
        fd.close()
        cmd="iterate d %s %s %s %s %s <  x.tmp" % (
            len(self.divs), K, maxdeg, n, mindeg)
        #cmd = os.path.join(os.path.split(__file__)[0], cmd)

        fd=os.popen(cmd)
        fd_arr=fd.readlines()
        fd.close()

        return(numpy.array(map(float, string.split(fd_arr[0]))))

    def degdist_numpy(self, maxdeg, n, mindeg=0):
        """Returns the degree distribution from mindeg to maxdeg degree.

        It should be used with the (iterated) link probability measure.

        Parameters:
            maxdeg: integer
                the maximal degree for we calculate the degree distribution
            n: integer
                the size of the network (the number of vertices)
            mindeg: integer
                the minimal degree for we calculate the degree distribution

        Returns:
            rho: the degree distribution as a list with length of maxdeg+1.
                The index k gives the probability of having degree k.

        Below mindeg the values of the return value will be zeros::

            >>> pm = ProbMeasure(m=3)
            >>> lpm = pm.iterate(K=2)
            >>> lpm.degdist(maxdeg=8, n=2000)
            array([  1.89094745e-11,   5.06324548e-10,   6.00749549e-09,
                     4.87734537e-08,   2.99004276e-07,   1.47045519e-06,
                     6.03452582e-06,   2.12437135e-05,   6.54696931e-05])
            >>> lpm.degdist(maxdeg=8, n=2000, mindeg=3)
            array([  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
                     4.87734537e-08,   2.99004276e-07,   1.47045519e-06,
                     6.03452582e-06,   2.12437135e-05,   6.54696931e-05])
        """

        assert isinstance(n, int)
        assert isinstance(maxdeg, int)
        assert n > maxdeg
        divs = self.divs
        n_intervals = len(divs)
        lengths = [divs[i] - divs[i-1] for i in xrange(1, n_intervals)]
        lengths.insert(0, divs[0])
        # Eq. 5, where ...
        avgdeg = [n*sum([self.probs[i][j]*lengths[j]
                      for j in xrange(n_intervals)])
                      for i in xrange(n_intervals)]
        log_avgdeg = numpy.log(avgdeg)
        rho = numpy.zeros(maxdeg+1)
        log_lengths = numpy.log(lengths)
        const = 0.5*math.log(2*math.pi)
        log_factorial = [const + (k+.5)*math.log(k) - k
                         for k in xrange(1, maxdeg+1) ]
        log_factorial.insert(0, 0)
        # Eq. 4
        for i in xrange(n_intervals):
            # Eq. 5
            log_rho_i = [(k * log_avgdeg[i] - log_factorial[k]  -  avgdeg[i])
                         for k in xrange(mindeg, maxdeg+1)]
            log_rho_i_length = numpy.array(log_rho_i) + log_lengths[i]
            #TODO the [mindeg:] slice is good? Should return mindeg as well?
            rho[mindeg:] += numpy.exp(log_rho_i_length)
        return rho

    degdist = degdist_numpy

    def degdist_bigfloat(self, maxdeg, n, mindeg=0):
        """Returns the degree distribution from 0 to maxdeg degree.

        It should be use with the (iterated) link probability measure.

        Parameters:
          maxdeg: the maximal degree for we calculate the degree distribution
          n: the size of the network (the number of vertices)

        Returns:
            rho: the degree distribution as a list with length of maxdeg+1.
                The index k gives the probability of having degree k.
        """

        assert isinstance(n, int) and isinstance(maxdeg, int) and n > maxdeg
        import bigfloat
        context = bigfloat.Context(precision=10)
        divs = self.divs
        n_intervals = len(divs)
        lengths = [divs[i] - divs[i-1] for i in xrange(1, n_intervals)]
        lengths.insert(0, divs[0])
        log_lengths = numpy.log(lengths)
        # Eq. 5, where ...
        avgdeg = [bigfloat.BigFloat(n*sum([self.probs[i][j]*lengths[j]
                  for j in xrange(n_intervals)]), context=context) for i in xrange(n_intervals)]
        #log_factorial = [ 0.5*bigfloat.log(2*math.pi, context=context) + (d+.5)*bigfloat.log(d, context=context) - d
        #                 for d in xrange(1,maxdeg+1) ]
        log_factorial = [bigfloat.log(bigfloat.factorial(k), context=context)
                         for k in xrange(1, maxdeg+1)]
        log_factorial.insert(0, 0)

        rho = [bigfloat.BigFloat(0, context=context)] * (maxdeg+1)
        # Eq. 4
        for i in xrange(n_intervals):
            # Eq. 5
            log_rho_i = [(bigfloat.mul(k, bigfloat.log(avgdeg[i]), context=context) - log_factorial[k] - avgdeg[i])
                         for k in xrange(mindeg, maxdeg+1)]
            log_rho_i_length = [log_rho_i[k] + log_lengths[i]
                         for k in xrange(mindeg, maxdeg+1)]
            for k in xrange(mindeg, maxdeg+1):
                rho[k] += bigfloat.exp(log_rho_i_length[k], context=context)
        return rho

    def degdist_mult(self, maxdeg, n):
        """Returns the degree distribution from 0 to maxdeg degree.

        Version with subsequent multiplying of the k=0 value.
        It will not work for n>1000 properly.

        It should be use with the (iterated) link probability measure.

        Parameters:
          maxdeg: the maximal degree for we calculate the degree distribution
          n: the size of the network (the number of vertices)

        Returns:
          rho: the degree distribution as a list with length of maxdeg+1.
            The index k gives the probability of having degree k.
        """

        assert isinstance(n, int) and isinstance(maxdeg, int) and n > maxdeg
        divs = self.divs
        n_intervals = len(divs)
        lengths = [divs[i] - divs[i-1] for i in xrange(1, n_intervals)]
        lengths.insert(0, divs[0])
        # Eq. 5, where ...
        avgdeg = [sum([self.probs[i][j]*lengths[j]
                    for i in xrange(n_intervals)]) for j in xrange(n_intervals)]
        avgdeg = [n*a for a in avgdeg]
        rho = []
        rho_i = [math.exp(-avgdeg[i]) for i in xrange(n_intervals)]
        for k in xrange(1, maxdeg+1):
            # Eq. 4
            rho.append( sum(rho_i[i] * lengths[i] for i in xrange(n_intervals)) )
            # Eq. 5
            rho_i = [(rho_i[i] * avgdeg[i] / k) for i in xrange(n_intervals)]
        rho.append( sum(rho_i[i] * lengths[i] for i in xrange(n_intervals)) )
        return rho

    def generate(self, n=200):
        """Creates a network according to the probability matrix.

        n: positive integer
           The number of vertices in the network.

        """
        network = igraph.Graph(n)
        values = numpy.random.random(n)
        parts = numpy.array([bisect_left(self.divs, val)
                   for val in values])
        edges = []
        for i, part_i in enumerate(parts):
            probs = self.probs[part_i, parts[i+1:]]
            links = numpy.random.random(len(probs)) < probs
            edges.extend((i, jj+i+1) for jj, linked in enumerate(links) if linked)
        network.add_edges(edges)
        return network

    def __iterate_divs(self, K=3):
        """Create iterated divs for self.iterate."""
        assert isinstance(K, int) and K > 0
        if K == 1:
            return self.divs
        olddivs = numpy.concatenate(([0], self.__iterate_divs(K-1)))
        newdivs = []
        for i in xrange(len(olddivs)-1):
            start, stop = olddivs[i], olddivs[i+1]
            interval_length = stop - start
            for k in self.divs:
                newdivs.append(start + interval_length*k)
        return newdivs

    def __iterate_probs(self, K=3):
        """Create iterated probabilities for self.iterate."""
        assert isinstance(K, int) and K > 0
        if K == 1:
            return self.probs
        m = self.m
        indices = numpy.arange(m**K)
        indices_dict = {}
        for i in indices:
            indices_dict[i] = []
        for iteration in xrange(K):
            indices2 = (indices % (m**(iteration+1))) // (m**iteration)
            for i in indices:
                indices_dict[i].append(indices2[i])

        probs = numpy.ones((m**K, m**K))
        for i in indices:
            ilist = indices_dict[i]
            for j in indices:
                jlist = indices_dict[j]
                for it in xrange(K):
                    probs[i,j] *= self.probs[ilist[it], jlist[it]]
        return probs

    def iterate(self, K=3, maxdivs=250):
        """Create a multifractal network with (K-1) iteration."""
        m = len(self.divs)
        divnum = m**K
        if divnum > maxdivs:
            # FIXME: other exception
            raise ValueError("There would be too many division points in the iterated network.\n"
                  "Division points in the wanted iteration {divnum}, the limit maxdivs={maxdivs}.\n"
                  "If you know, what do you do, set the maxdivs parameter of iterate method higher."
                  .format(divnum=divnum, maxdivs=maxdivs))
        divs, probs = self.__iterate_divs(K), self.__iterate_probs(K)
        return ProbMeasure(divs, probs, output=self.out, iterated=True)

    def newdivs(self, n=2, ii=None):
        """Change one of the divs.

        Parameters:
            n: integer, default 2
              The exponent in the calculation of the new value of the divs.
            ii: integer in [0, m-1] or None, default None
              The index of the divs.
              If None, it chooses randomly.


        New value of the d[i] is:

            /     a \ n
            |p - ---|
            \     l /                      a
        l * ---------     + d[i],  if p > ---
            /     a \ n-1                  l
            |1 - ---|
            \     l /

            /     a \ n
            |p - ---|
            \     l /                      a
        l * ---------     + d[i],  if p < --
            /     a \ n-1                  l
            |  - ---|
            \     l /

        where "p" is a random variate from a uniform distribution in [0,1[
        "l" and "a" are drawn below:

        |<----------- l ----------->|
        |<---- a ----->|
        |--------------|------------|
        d[i-1]        d[i]          d[i+1]

        """
        assert isinstance(n, int) and n > 0, "the n exponent must be positive integer"
        if ii is None:
            ii = randrange(self.m-1)
        divs = self.divs[:]
        if ii == 0:
            l = divs[ii+1]
            a = divs[ii]
        else:
            l = divs[ii+1] - divs[ii-1]
            a = divs[ii]   - divs[ii-1]
        p = random()
        if p > a/l:
            delta = l*(p-a/l)**n/(1-a/l)**(n-1)
        else:
            delta = l*(p-a/l)**n/( -a/l)**(n-1)
        #TODO with delta = -l*(-p+a/l)**n / abs(a/l)**(n-1) we can use rational exponents
        # We need to test, whether it is equivalent
        # Why we need abs?
        #TODO a separate delta function to make unittests easier.
        divs[ii] += delta
        assert 0 < divs[ii] < 1, "divs[{0}] = {1} not in ]0, 1[."
        self.out.write("divs changed to %s" % divs)
        return divs

    def newprobs(self, maxrel=0.1):
        """Change one of the probabilities.

        And change the others too, to make the sum 1.

        Parameters:
            maxrel: float in ]0,1], default 0.1 (10  %)
                maximal relative change
        """
        m = self.m
        ii = randrange(m)
        jj = randrange(m)
        self.out.write("i=%d, j=%d" % (ii, jj))
        probs = self.probs.copy()
        p = random()
        delta = probs[ii, jj] * (2 * p - 1) * maxrel
        probs[ii, jj] += delta
        if ii != jj:
            probs[jj, ii] += delta
        probs /= numpy.sum(probs)
        assert (numpy.sum(probs) - 1) < 1e-6
        self.out.write("probs changed to %s" % probs)
        return probs

#from cyProbMeasure import cyProbMeasure as ProbMeasure

def simpleProbMeasure(m, **kwargs):
    """Creates a ProbMeasure with equal probabilities and equal interval lengths on axes.

    Parameters:
        m: integer
            the number of intervals on each axes
    """
    assert m > 1 and isinstance(m, int)
    divs = numpy.linspace(1, 0, m, endpoint=False)
    divs = divs[::-1]
    probs = 1. / m * numpy.ones((m, m))
    return ProbMeasure(divs, probs, **kwargs)

def get_Tfactor(q, steps):
    """Gives back Tfactor for a (Tlimit/T0, steps) pair"""
    assert steps % 2 == 0, "steps is always even"
    return (q*0.9999999)**(1/(steps/2))

def steps(T0, Tfactor, Tlimit):
    """Returns with the number of steps

    Parameters (see at the Generator class):
        T0:
            initial temperature
        Tfactor:
            temperature factor
        Tlimit:
            temperature limit

    """
    return 2 * int(numpy.ceil(numpy.log(Tlimit/T0) / numpy.log(Tfactor)))

class Generator(object):
    """Generate a network with given properties.

    Paremeters:
        T: float
            initial temperature
        steps: int
            the number of steps
        Tlimit: float
            the temperature when the generation stops
        m: integer
            the probmeasure will be mxm type
        K: integer
            it will use K iteration to average the energy
        divexponent: int
            the exponent in the formula for adjusting division points
        project: string
            The values will be stored in the directory named with
            ``'project_' + project``
            E.g if project is ``'base'`` the directory will be ``'project_base'``.

    The generation of a generator measure - an example::

        generator = mfng.Generator(T0=0.2, steps=10000, Tlimit=0.00002,
                m=3, K=3,
                n=8000,
                divexponent = 7,
                project = 'base',
                bigfloat = False
                division = True,
                )
        generator.append_property(
                mfng.AverageDegree(avg_degree)
                )
        generator.go()
    """
    def __init__(self,
            T0, steps, Tlimit,
            m=2, K=3,
            n=200,
            divexponent=2,
            project = "base",
            verbose = False,
            ):
        """Set parameters and make the initial probmeasure.
        """
        self.project = project
        self.project_dir = "project_%s" % project
        if not os.path.isdir(self.project_dir):
            os.mkdir(self.project_dir)
        #self.out = Output(os.path.join(self.project_dir, "details.txt"), with_stdout=verbose) #TODO
        self.out = NoOutput() #os.path.join(self.project_dir, "details.txt"), with_stdout=verbose)
        self.T = self.T0 = T0
        self.Tfactor = get_Tfactor(Tlimit/T0, steps)
        print("Tfactor =", self.Tfactor)
        assert 0 < self.Tfactor < 1
        self.Tlimit = Tlimit
        self.K = K
        self.n = n
        assert isinstance(divexponent, int) and divexponent > 0, "divexponent must be positive integer"
        self.divexponent = divexponent
        self.properties = []
        self.accept_reject = ""

        # Should be greater then the energy of the initial probmeasure.
        self.E = self.Ei = 1e310
        assert isinstance(m, int) and m > 1
        self.probmeasure = simpleProbMeasure(m, output=self.out)

        self.starttime = self.stoptime = ""
        self.verbose = verbose # If True prints energies and temperatures.
        self.out.write("=============================================\n"
              "The properties of the Generator object"
              )
        self.out.write(self.settings())
        self.out.write("date = %s" % hms(time()))
        if not isinstance(steps, int):
            steps = int(steps)
        if steps % 2 != 0:
            steps += 1
            print("The number of steps must be even. The generation will have {0} steps".format(steps))
        self._steps = steps
        print("\nThere will be %d steps." % self._steps)
        self.energy_list = []
        self.divs_list = []
        self.probs_list = []

    def append_property(self, prop):
        """Appends the property prop to the list of properties."""
        self.properties.append(prop)
        assert prop.probmeasure_type in ["iterated", "generator"],\
                "probmeasure_type should be 'iterated' or 'generator'"
        self.out.write("\n%s added." % prop)

    def has_accepted(self, probmeasure):
        """Averages the energy of generated networks.

        If the probmeasure has accepted, it saves the new energy and the
        accepted probmeasure.

        Paramters:
            probMeasure:
                the base ProbMeasure instance

        Returns
            True if the new probmeasure has accepted, else not.

        """

        energies = []
        values = []
        E = 0.0
        lpm = False
        for prop in self.properties:
            if prop.probmeasure_type == "iterated":
                if lpm is False:
                    lpm = probmeasure.iterate(self.K)
                probm = lpm
            else:
                probm = probmeasure
            E_, value = prop.energy(probm, self.n)
            values.append((prop.shortname, value))
            E += E_
        self.out.write(" %s\n" % ",".join("{0}={1}".format(s, v) for s, v in values))
        has_accepted = False
        result = None
        Eorig = self.E
        if E < self.E:
            has_accepted = True
            result = "+++ ACCEPTED +++"
            self.accept_reject += "A"
            print("A") #TODO
        elif random() < numpy.exp(-(E-self.E)/self.T):
            has_accepted = True
            result = "+++ ACCEPTED in the E_new > E_original branch. +++"
            self.accept_reject += "a"
            print("a") #TODO
        if has_accepted:
            self.E = E
            self.probmeasure = probmeasure
        if not result:
            result = "--- REJECTED ---"
            self.accept_reject += "."
            print(".") #TODO
        self.out.write("  E=%10s (orig. %7s)   %s\n" % (E, Eorig, result))
        return has_accepted

    def doublestep(self):
        """One double step includes one newprobs and one newdivs method.
        Temperature is adjusted in each double step.
        """
        self.out.write("T=%s" % self.T)

        for change in ["probs", "divs"]:
            self.step += 1
            pm = self.probmeasure
            if change == "probs":
                self.out.write(" Changing the probabilities. ".center(42, "*"))
                probs = pm.newprobs()
                divs = pm.divs
            else:
                self.out.write("  Changing the divs.  ".center(42, "*"))
                probs = pm.probs
                divs = pm.newdivs(n=self.divexponent)

            self.energy_list.append(self.E)
            self.divs_list.append(self.probmeasure.divs)
            self.probs_list.append(self.probmeasure.probs)

            if self.probmeasure.m == 3:
                print("divs=[{0[0]:6.4f} {0[1]:6.4f}] p00={probs[0][0]:6.4f} "
                "p01={probs[0][1]:6.4f} p02={probs[0][2]:6.4f} p11={probs[1][1]:6.4f} "
                "p12={probs[1][2]:6.4f} p22={probs[2][2]:6.4f} "
                "T={T:3e} E={E:.8f}".format(divs, probs=probs, T=self.T, E=self.E), end=" ") #TODO
            elif self.probmeasure.m == 2:
                print("div={0[0]:6.4f} p00={probs[0][0]:6.4f} "
                "p01={probs[0][1]:6.4f} p11={probs[1][1]:6.4f} "
                "T={T:3e} E={E:.8f}".format(divs, probs=probs, T=self.T, E=self.E), end=" ") #TODO
            else:
                print("divs={0:60} "
                "T={T:3e} E={E:.8f}".format(divs, probs=probs, T=self.T, E=self.E), end=" ") #TODO

            newpm = ProbMeasure(divs, probs, output=self.out)
            has_accepted = self.has_accepted(newpm)

        self.T *= self.Tfactor
        self.accept_reject += " "

    def go(self):
        """Starts the generation."""
        self.starttime = time()
        assert self.has_accepted(self.probmeasure)
        self.accept_reject = ""
        self.Ei = self.E
        self.step = 0
        while self.T > self.Tlimit:
            self.out.write("=============\nSteps %d & %d\n=============" % (self.step+1, self.step+2))
            self.doublestep()
        self.stoptime = time()
        self.out.write(self.results(with_accept=True))
        if self._steps != self.step:
            print("The calculated steps ({0}) not equals"
                  " with the occurred one ({0})."
                  "write to the author, please.".format(self._steps, self.step))
        self.write_run()

        return self.probmeasure

    def settings(self, compact=False):
        """Gives the settings as string"""

        if compact:
            form = (
            "     T0=%s,  Tfactor=%s, Tlimit=%s,\n"
            "     m=%s, K=%s,\n"
            "     n=%s,\n"
            "     divexponent=%s,\n"
            "     properties= [\n        %s\n     ],"
            )
        else:
            form = (
            "======================================\n"
            "T0 = %s,  Tfactor = %s, Tlimit = %s\n"
            "m=%s, K=%s\n"
            "n=%s\n"
            "divexponent=%s\n"
            "Properties:\n        %s\n"
            )
        pm = self.probmeasure

        properties = ",\n        ".join(['"%s"' % prop.name for prop in self.properties])
        return form % (
                self.T0, self.Tfactor, self.Tlimit,
                pm.m, self.K,
                self.n,
                self.divexponent,
                properties,
                )

    def results(self, with_settings=True, with_accept=False, compact=False):
        """Gives the result as string"""

        pm = self.probmeasure
        runtime = (self.stoptime - self.starttime)/60
        if compact:
            acceptform = "\n     accept_reject = \"%s\","
            form = (
                "     divs= %s,\n"
                "     probs=\\\n       %s,\n"
                "     #Energies\n       Ei =%10f, #Initial\n       Ef =%10f, #Final\n"
                "     runtime =  %f, # minutes (%s -- %s)\n"
                "     steps =  %d,\n"
                )
        else:
            acceptform = "Accept with Enew < E (A) or with Enew > E (a); reject (.):\n%s\n"
            form= (
                "______________________________________\n"
                "Division points and the probability matrix\n"
                "divs= %s\n"
                "probs=\\\n%s\n"
                "Energies\n Initial E=%s\n Final   E=%s\n"
                "Runtime %f min (%s -- %s)\n"
                "Steps =  %d\n"
                )
        results =  form % (
                pm.divs,
                repr(pm.probs),
                self.Ei,
                self.E,
                runtime, hms(self.starttime), hms(self.stoptime),
                self.step,
                )
        if with_settings:
            results = self.settings(compact=compact) + results
        if with_accept:
            results += acceptform % self.accept_reject
        return results


    def write_run(self):
        """Write settings and results as Python dictionary."""
        runsfile = os.path.join(self.project_dir, "runs.py")
        runsjson = os.path.join(self.project_dir, "runs.json")
        runsshelf = os.path.join(self.project_dir, "runs.shelf")
        is_runsfile = os.path.isfile(runsfile)
        is_runsshelf = os.path.isfile(runsshelf)
        if is_runsshelf:
            db = shelve.open(runsshelf)
            self.labels = db.keys()
            db.close()
        else:
            self.labels = []
        reserved_indexes = [int(l.split("_")[1]) for l in self.labels if l.startswith("%d_" % self.n)]
        if reserved_indexes:
            index =  max(reserved_indexes)+1
        else:
            index = 1
        label = "%d_%03d" % (self.n, index)
        self.labels.append(label)

        text = (
            "runs[\"%s\"] = dict(\n"
            "  # Settings\n"
            "%s\n"
            "  # Results\n"
            "%s\n\n"
            ")\n\n"
            ) % (
            label,
            self.settings(compact=True),
            self.results(compact=True, with_settings=False, with_accept=True)
            )

        rf = open(runsfile, "a")
        if not is_runsfile:
            rf.write("""#!/usr/bin/env python
# -*- coding: utf-8 -*-

try:
    from numpy import array
except ImportError:
    array = lambda x: x

runs={}
""")
        rf.write(text)
        rf.close()

        pm = self.probmeasure
        values_to_store = dict(
                # Settings
                T0=self.T0,
                Tfactor=self.Tfactor,
                Tlimit=self.Tlimit,
                m=pm.m,
                K=self.K,
                n=self.n,
                divexponent=self.divexponent,
                properties = [prop.name for prop in self.properties],
                # Results
                divs = pm.divs,
                probs = [list(probsrow) for probsrow in pm.probs],
                Ei = float(self.Ei),
                E = float(self.E),
                runtime = (self.stoptime - self.starttime)/60,
                starttime = self.starttime,
                stoptime = self.stoptime,
                steps = self.step,
                accept_reject = self.accept_reject,
                energy_list = self.energy_list,
                divs_list = self.divs_list,
                probs_list = self.probs_list,
                )
        key_order = """T0 Tfactor Tlimit m K n divexponent properties divs probs
                    Ei E runtime starttime stoptime steps accept_reject
                    energy_list probs_list divs_list
                    """.split()
        keys1, keys2 =  set(key_order), set(values_to_store.keys())
        if not  keys1 == set(values_to_store.keys()):
            print("Keys in key_order and values_to_store are not the same.")
            print("These values differs:",
                keys1.symmetric_difference(keys2))
        db = shelve.open(runsshelf)
        assert label not in db.keys()
        db[label] = values_to_store
        db.close()
        del values_to_store["probs_list"] # TODO should convert to list of lists
        json_text = json.dumps(values_to_store, indent=4)
        f = open(runsjson, "a")
        f.write('"%s" :\n' %label)
        f.write(json_text)
        f.write(',\n')
        f.close()

        self.out.write(
                "I have written the results with the label '%s' into the files:\n"
                " - '%s' (text format, Python code),\n"
                " - '%s' (binary format, Python shelf).\n"
                   % (label, runsfile, runsshelf,))

class DistributionFunction(object):
    """The main class to give a degree distribution as a function to the Generator.

    Parameters:
        funtion: string
          A valid Python expression with the only variable 'k'.
          This describes the target degree distribution.
          The function will be normalized.
          numpy functions are allowed.
          E.g.:   "5*numpy.exp(-4*k)"
        maxdeg: integer
          the maximal degree we compare it the degree distribution
          of the probability measure with.
        mindeg: int, default 0, 0 <= mindeg < maxdeg
          like maxdeg, but for minimal degree
    """
    def __init__(self, function, maxdeg, mindeg=0, **kwargs):
        assert isinstance(maxdeg, int)
        assert isinstance(mindeg, int)
        assert 0 <= mindeg < maxdeg
        self.probmeasure_type = "iterated"
        function = re.sub(r"\bd\b", "k", function)
        self.function_str = function
        exec("self.function = lambda k: {0}".format(function))
        self.shortname = "degdistfunc"
        self.name = "DistributionFunction('{0}', maxdeg={1}, mindeg={2})".format(
                self.function_str, maxdeg, mindeg)
        self.mindeg, self.maxdeg = mindeg, maxdeg
        # the values for all of the degrees in the interval [mindeg, maxdeg]
        self.degdist = numpy.array([0]*mindeg + [self.function(k) for k in xrange(mindeg, maxdeg*100+1)])
        assert min(self.degdist[mindeg:maxdeg+1]) > 0
        # normalisation
        self.degdist /= float(sum(self.degdist))
        self.degdist = self.degdist[:maxdeg+1]
        print("Sum degdist: {0}".format(sum(self.degdist)))#TODO
        self.division = division

        keys = set(kwargs.keys())
        if isinstance(self, DistributionFunctionC) and "K" in keys:
            keys -= set('K')
        if keys:
            print("These arguments are deprecated or erroneous: ", ", ".join(keys))

    def __str__(self):
        return self.name

    def __repr__(self):
        return self.name

    def value(self, k):
        return self.function(k)

    def energy(self, link_probmeasure, n):
        """Returns the energy of the link probability measure.

        Parameters:
            link_probmeasure: ProbMeasure
                The probability measure we want to measure the energy for.
            n: int, n > maxdeg
                The size of the networks (number of vertices), we calculate the
                degree distribution from the link probability measure for.

        """
        assert isinstance(n, int) and n > self.maxdeg
        # Eq. 6
        actual_degdist = link_probmeasure.degdist_numpy(self.maxdeg, n)
        sumdd = sum(actual_degdist)
        degdist = self.degdist
        #summa = sum(abs(actual_degdist[k]/degdist[k] - 1)
        #            for k in range(self.mindeg, self.maxdeg+1))
        summa = sum(abs(degdist[k] - actual_degdist[k])/max(degdist[k], actual_degdist[k])
                for k in xrange(self.mindeg, self.maxdeg+1))
        return float(summa), sumdd

class DistributionFunctionC(DistributionFunction):
    """Faster version of DistributionFunction written in C++

    This program iterates the ProbMeasure itself.
    """
    def __init__(self, *args, **kwargs):
        super(DistributionFunctionC, self).__init__(*args, **kwargs)
        self.probmeasure_type = "generator"
        self.K = kwargs.get("K")

    def energy(self, probmeasure, n):
        """Returns the energy of the link probability measure.

        Parameters:
            probmeasure: ProbMeasure
                The probability measure we want to measure the energy for.
            n: int, n > maxdeg
                The size of the networks (number of vertices), we calculate the
                degree distribution from the link probability measure for.

        """
        assert isinstance(n, int) and n > self.maxdeg
        # Eq. 6
        actual_degdist = probmeasure.degdist_iterated(self.maxdeg, n)
        sumdd = sum(actual_degdist)
        degdist = self.degdist
        summa = sum(abs(actual_degdist[k]/degdist[k] - 1)
                    for k in range(self.mindeg, self.maxdeg+1))
        #summa = sum(abs(degdist[k] - actual_degdist[k])/max(degdist[k], actual_degdist[k])
        #        for k in xrange(self.mindeg, self.maxdeg+1))
        return float(summa), sumdd

def log_binner(num_list):
    """Returns the log-binned list.

    The length of a bin is two times the previous bin.
    It does not divide with the length of the bin, and does not check,
    whether the sum of the values is 1 (or less then 1).

    Example::

        >>> log_binner([0.5, 0.25, 0.125, 0.125])  #doctest: +SKIP
        [0.5, 0.375, 0.125]
        >>> log_binner([1]*7)                      #doctest: +SKIP
        [1, 2, 4]

    """
    if not isinstance(num_list, (list, tuple, numpy.ndarray)):
        raise ValueError("The argument of log_binner must be list or tuple.")
    for num in num_list:
        if not isinstance(num, (int, long, float)):
            raise ValueError("The elements of the argument of log_binner must be numbers.")
    binned_list = []
    i = 1
    max_index = len(num_list)
    while i <= max_index:
        if i*2 <= max_index:
            binned_list.append(sum(num_list[i-1:i*2-1]))
        else:
            binned_list.append(sum(num_list[i-1:]))
        i *= 2
    return binned_list

class LogBinnedDegDist(object):
    """The class to give a logarithmic binned degree distribution to the Generator.

    Its main role to create a network with a degree distribution like that of
    an original network.

    Parameters:
        degrees: list
          The list of degrees of a network.
    """
    def __init__(self, degrees):
        assert isinstance(degrees, list)
        self.probmeasure_type = "iterated"
        for deg in degrees:
            assert isinstance(deg, (int, long)) and deg >= 0
        self.shortname = "LogBinnedDegDist"
        self.name = "LogBinnedDegDist('{0}')".format(degrees)
        self.n = n = len(degrees)
        k_max = max(degrees)
        degdist = [degrees.count(k) for k in xrange(1, k_max+1)]
        binned_degdist = log_binner(degdist)
        self.binned_degdist = numpy.array(binned_degdist)/n
        self.maxdeg = 2**len(binned_degdist) - 1

    def __str__(self):
        return self.name

    def __repr__(self):
        return self.name

    def value(self):
        return self.binned_degdist

    def energy(self, link_probmeasure, n):
        actual_binned_degdist = log_binner(link_probmeasure.degdist(self.maxdeg, n, mindeg=1)[1:])
        actual_binned_degdist = numpy.array(actual_binned_degdist)
        energy = 0
        sum_bdd = sum(actual_binned_degdist)
        diffsum = sum(abs(self.binned_degdist - actual_binned_degdist)/ \
                      numpy.maximum(self.binned_degdist,  actual_binned_degdist)
                )
        return 1.*diffsum, sum_bdd

# TODO  Not ready.
class RealLogBinnedDegDist(object):
    def __init__(self, archive_file):
        import cxnet
        archive, file_ = cxnet.archives.get_archive_name(archive_file)
        file_name = os.path.join("..", "netdata", "{0}.gml".format(file_))


class AverageDegree(object):
    """The class to give an average degree to the Generator.

    Parameter:
        value: float, value >= 0
          The target average degree.
    """
    def __init__(self, value):
        assert value >= 0
        self.probmeasure_type = "iterated"
        self._value = value
        self.name = "AverageDegree({0})".format(value)
        self.shortname = "avgdeg"

    def __str__(self):
        return self.name

    def __repr__(self):
        return self.name

    def value(self):
        return self._value

    def energy(self, link_probmeasure, n):
        """Returns the energy of the link probability measure.

        Parameter:
            link_probmeasure: ProbMeasure
                The probability measure we want to measure the energy for.
            n: int
                The size of the networks (number of vertices), we calculate the
                degree distribution from the link probability measure for.

        """
        assert isinstance(link_probmeasure, ProbMeasure) and isinstance(n, int)
        actual_avg_degree = link_probmeasure.avg_degree(n)
        energy = abs(self._value - actual_avg_degree)/max(self._value, actual_avg_degree)
        return energy, self._value

if __name__ == '__main__':
    import doctest
    doctest.testmod()
