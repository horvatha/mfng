#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Classes and function to analyze data of the results of the mfng projects.

The Probmeasure in this file is the derived class that of the mfng class.
It uses pylab to plot the measure.
"""

import mfng
from mfng import Output

try:
    import pylab
except RuntimeError:
    import matplotlib
    matplotlib.use("Agg")
    import pylab
import cxnet
import sys
import os
import shelve
import json
sys.path.insert(0, ".")

class NoProjectError(Exception):
    def __init__(self, project):
        self.project = project
    def __str__(self):
        return "There is no project named '%s'." % self.project

class EmptyProjectError(Exception):
    def __str__(self):
        return "There is no generated measure in this project."

class NoLabelsError(Exception):
    def __str__(self):
        return "There is no labels. Set it with the set_labels() method."

class ProbMeasure(mfng.ProbMeasure):
    """The Probmeasure in this file is the derived class that of the mfng class.

    It uses pylab to plot the measure.
    """
    def plot(self, ticklabels=True, colorbar=True, **kwargs):
        """Plots the Probability measure.

        Parameters:
        - `ticklabels`: if True [default] it prints the values at the axes
        - `kwargs`: keyword argumentums to pass to the pcolor function of pylab

        """
        X=Y=pylab.array([0] + self.divs)
        fig = pylab.figure()
        ax = fig.add_subplot(111)
        #if "cmap" not in kwargs:
        #    kwargs["cmap"] = pylab.cm.Blues
        pylab.pcolor(X,Y,self.probs, **kwargs)
        pylab.axis((0,1,1,0))
        title = "Link probability measure" if self.iterated else "Generating measure"
        pylab.suptitle(title)
        # Perhaps in matplotlib >1.1 should use title instead of suptitle and
        # matplotlib.pyplot.tight_layout
        if colorbar:
            pylab.colorbar()
        for tick in ax.xaxis.get_major_ticks():
            if ticklabels:
                tick.label2On=True
            tick.label1On=False
        if not ticklabels:
            for tick in ax.yaxis.get_major_ticks():
                tick.label1On=False
        return ax

    def iterate(self, K=3, maxdivs=730):
        """Create a multifractal network with (K-1) iteration.

        We need to redefined to give back the Probmeasure module with plot in this modul.
        """
        divnum = len(self.divs)**K
        if divnum > maxdivs:
            raise ValueError("There would be too many division points in the iterated network.\n"
                  "Division points in the wanted iteration {divnum}, the limit maxdivs={maxdivs}.\n"
                  "If you know, what do you do, set the maxdivs parameter of iterate method higher."
                  .format(divnum=divnum, maxdivs=maxdivs))
        divs, probs = self.__iterate_divs(K), self.__iterate_probs(K)
        return ProbMeasure(divs, probs, output=self.out, iterated=True)


def avg_sigma(numlist):
    """Returns with the average and sigma of the listed numbers.

    Parameter:
      `numlist`: list or array

    """
    average = pylab.average(numlist)
    numlist = pylab.array(numlist)
    sum2 = sum((numlist - average)**2)
    sigma = pylab.sqrt(sum2 / (len(numlist)-1))
    return average, sigma


class Runs:
    """Class for analyzing the runs.

    Parameter:
      `project`: None [default] or the name of the project.
          If None it will search the ``base`` project. If it does not find,
          it will use the actual directory as the project directory.

    """

    def __init__(self, project=None):
        if project is None:
            if os.path.isdir("project_base"):
                project_dir = "project_base"
            else:
                project_dir = "."
        elif isinstance(project, str) and os.path.isdir(project):
            project_dir = project
        elif not isinstance(project, str):
            raise ValueError("project must be string or None.")
        else:
            project_dir = "project_%s" % project
        if not os.path.isdir(project_dir):
            raise NoProjectError(project)
        self.project_dir = project_dir
        if isinstance(project, dict):
            runs = project
        else:
            db = shelve.open(os.path.join(project_dir, "runs.shelf"))
            runs = {}
            for key in db:
                runs[key] = db[key]
            db.close()
        self.runs = runs
        self.labels = None
        if not self.runs.keys():
            raise EmptyProjectError

    def load_from_json(self):
        """Experimental loader."""
        f = open(os.path.join(self.project_dir, "runs.json"))
        text = "".join(f.readlines())
        f.close()
        text = "{%s\n}" % text[:-2]
        self.runs2 = json.loads(text)

    def set_labels(self, labels="all"):
        """Set the labels for the analyzis.

        If we add the list of labels, it keeps the order.
        Otherwise the labels will be sorted.

        Parameter:

            `labels`: string or list of strings
               This can be:

               * 'all' will include all labels available,
               * one label we want to analyze,
               * the start of the labels we want to analyze,
               * the list of the labels we want to analyze.

        Returns:

            With the labels it have set.

        Example::

            >>> sorted(r.runs.keys())  # All labels
            ['2000_007', '2000_008', '2000_009', '2000_010', '2000_011']
            >>> r.set_labels("all")
            ['2000_007', '2000_008', '2000_009', '2000_010', '2000_011']
            >>> r.set_labels("2000_008")
            ['2000_008']
            >>> r.set_labels("2000_01")
            ['2000_010', '2000_011']
            >>> r.set_labels(["2000_009", "2000_007"])
            ['2000_009', '2000_007']

        """
        all_labels = self.runs.keys()
        if not isinstance(labels, list):
            if labels == "all":
                labels = all_labels
            elif labels in all_labels:
                labels = [labels]
            elif isinstance(labels, str):
                labels = [ l for l in all_labels if l.startswith(labels) ]
            labels.sort()
        assert labels and all((label in all_labels) for label in labels )
        self.labels = labels
        return labels

    def properties(self, n=100):
        """
        Prints the statistical properties (maximal and average degree) of the results.

        It goes through the labels we have set with `self.set_labels`.

        Parameters:

            `n`: optional
                the number of the generated networks to create statistics.

        Example::

            >>> r.set_labels("200_020")
            ['200_020']
            >>> r.properties(n=100)
            ====================
            label = 200_020
            number of generated networks = 100
            divs = [0.81104014013270886, 1.0],
            probs=[[ 0.27609856  0.23149258]
             [ 0.23149258  0.26091629]]
               avg max deg =  9.87+-1.0115993937, avg average deg=3.5651+-0.175231419857

        """
        out = Output(os.path.join(self.project_dir, "properties.txt"))
        assert isinstance(n, int) and n>0
        if not isinstance(self.labels, list) or not self.labels:
            raise NoLabelsError
        for label in self.labels:
            run = self.runs[label]
            doubleline = "=" * 20
            out.write(doubleline)
            out.write("label = %s" % label)

            divs, probs = run["divs"], run["probs"]
            pm = mfng.ProbMeasure(divs, probs)
            K = run.get("K") or run.get("k")  # Older files have k instead of K.
            ipm = pm.iterate(K=K)
            out.write("number of generated networks = %s\ndivs = %s,\nprobs=%s" % (n,divs, probs))
            maxdegs = []
            avgdegs = []
            for i in range(n):
                nw = ipm.generate(n=run["n"])
                deg = nw.degree()
                maxdegs.append(max(deg))
                avgdegs.append(pylab.average(deg))
            avgmaxdeg, avgmaxdegsigma = avg_sigma(maxdegs)
            avgavgdeg, avgavgdegsigma = avg_sigma(avgdegs)
            out.write("\navg max deg = %5s+-%5s, avg average deg=%5s+-%5s\n" %
                    (avgmaxdeg, avgmaxdegsigma, avgavgdeg, avgavgdegsigma))
        del out

    def degdist(self, label=None, initial=False):
        """Returns with the degdist object."""
        if label is None:
            if self.labels is None:
                raise NoLabelsError
            label = self.labels[-1]
        assert label in self.runs.keys()
        ipm = self.probmeasure(initial=initial, iterated=True, label=label)
        run=self.runs[label]
        nw = ipm.generate(n=run["n"])
        dd = cxnet.DegreeDistribution(nw)
        return dd

    def loglog(self, with_initial=True, binning="ondemand"):
        """Plots the degree distribution of a generated network for each label.

        Parameters:
            `with_initial`: boolean, default True
              If it is True [default] it plots the initial degree distribution
              of the results with the first label in self.labels.
            `binning`: string
              The type of binning you want to use (ondemand, all or log).

        """
        if self.labels is None:
            raise NoLabelsError
        for label in self.labels:
            dd = self.degdist(label)
            dd.set_binning(binning)
            dd.loglog(label=label)
        if with_initial:
            ddi = self.degdist(self.labels[0], initial=True)
            ddi.set_binning(binning)
            ddi.loglog(label="initial")
        pylab.legend()

    def loglog1(self, with_initial=True, binning="ondemand", **kwargs):
        """Plots the degree distribution of a generated network for one label.

        This label come from the label argument or from the last item of the
        self.labels.
        """
        assert isinstance(with_initial, (bool, int)), "with_initial should be True, False, 1 or 0"
        label = kwargs.get("label")
        if label is None:
            if self.labels is None:
                raise NoLabelsError
            else:
                label=self.labels[-1]
        run = self.runs[label]
        if "K" not in run.keys():
            run["K"] = run["k"]  # Older files have k instead of K.
        n = run["n"]
        properties = run["properties"]
        if len(properties)==1 and properties[0].startswith("DistributionFunction"):
            from mfng import DistributionFunction
            exec("df={0}".format(properties[0]))
            pylab.plot(df.degdist, "b:", label="target dist.")
            maxdeg, mindeg = df.maxdeg, df.mindeg
            title= 'Deg. dist. {label} ' +\
                '(m={m},K={K},n={n},maxdeg={maxdeg},steps={steps})'.format(
                    maxdeg=maxdeg, label=label, **run)
        else:
            mindeg = 1
            maxdeg = int(0.03*n) # TODO What should be maxdeg?
            title= 'Deg. dist. {label} ' +\
                '(m={m},K={K},n={n},steps={steps})'.format(
                    label=label, **run)
        if "mindeg" in kwargs:
            mindeg = kwargs.get("mindeg")
        if "maxdeg" in kwargs:
            maxdeg = kwargs.get("maxdeg")

        dd = self.degdist(label)
        dd.set_binning(binning)
        dd.loglog(label="final gen.", markerfacecolor="g", markeredgecolor="g")

        ipm = self.probmeasure(iterated=True, initial=False, **kwargs)
        dd = ipm.degdist(maxdeg=maxdeg, n=n, mindeg=mindeg)
        pylab.plot(dd, "g--", label="final calc.")

        if with_initial:
            ddi = self.degdist(label, initial=True)
            ddi.set_binning(binning)
            ddi.loglog(label="initial gen.", markerfacecolor="red", marker="+")

            ipm = self.probmeasure(iterated=True, initial=True, **kwargs)
            dd = ipm.degdist(maxdeg=maxdeg, n=n, mindeg=mindeg)
            pylab.plot(dd, "r:", label="initial calc.")

        pylab.legend()
        axis_ = list(pylab.axis())
        axis_[2]=.01/n
        pylab.axis(axis_)
        pylab.title(title)
        return label

    def probmeasure(self, iterated=True, initial=False, **kwargs):
        """Returns the link probability or the generator measure.

        It can be plotted then.

        Parameters:
            `iterated`: boolean
                If True it returns with the iterated measure.
                (with the link probability measure)
            `initial`: boolean
                If True it returns with the initial measure.
            `label`: string or None
                If string, it will use this as label.

        Example. In this example the link probability measure and the generator
        measure will be plotted. ::

            >>> pm=r.probmeasure()
            >>> pm.plot(ticklabels=False)
            <matplotlib.axes.AxesSubplot object at 0x8f05d0c>
            >>> pm=r.probmeasure(iterated=False)
            >>> pm.plot()
            <matplotlib.axes.AxesSubplot object at 0x9632eec>

        """
        label = kwargs.get("label")
        if not label and (not isinstance(self.labels, list) or not self.labels):
            raise NoLabelsError
        if not label:
            label=self.labels[-1]
        run = self.runs[label]

        if initial:
            pm = mfng.simpleProbMeasure(m=run["m"])
        else:
            divs = run["divs"]
            probs = run["probs"]
            pm = ProbMeasure(divs, probs)

        if iterated:
            K = run.get("K") or run.get("k")  # Older files have k instead of K.
            pm = pm.iterate(K)

        return pm

    def plot_energy_list(self, **kwargs):
        """Plots the energies in the steps of the generation.

        It plots for all the labels in the self.labels.
        """
        if not isinstance(self.labels, list) or not self.labels:
            raise NoLabelsError
        number_of_plots = 0
        for label in self.labels:
            run = self.runs[label]
            if "energy_list" not in run:
                print("Sorry, there was no energy_list stored in the run '{0}'.".format(label))
                continue
            number_of_plots += 1
            pylab.plot(range(1,len(run["energy_list"])+1), run["energy_list"], label=label, **kwargs)
        if number_of_plots>1:
            pylab.legend(loc=1)
        pylab.title("Changing of energy")
        pylab.xlabel("step")
        pylab.ylabel("energy")

    def plot_divs_list(self, label=None, **kwargs):
        """Plots the divs in the steps of the generation.

        It plots for the first last in the self.labels.
        """
        if not label:
            if not isinstance(self.labels, list) or not self.labels:
                raise NoLabelsError
            else:
                label = self.labels[-1]
        run = self.runs[label]
        if "divs_list" not in run:
            print("Sorry, there was no divs_list stored in the run '{0}'.".format(label))
            return
        divs_list=run["divs_list"]
        for i in range(len(divs_list[0])):
            divs_i = [divs[i] for divs in divs_list]
            pylab.plot(range(1,len(divs_list)+1), divs_i, **kwargs)
        pylab.title("Changing of division points")
        pylab.xlabel("step")
        pylab.ylabel("division points")


if __name__ == "__main__":
    pass



