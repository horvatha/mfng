======================
cmfng Python package
======================

A realization of Multifractal Network Generator.

See:
Palla - Vicsek - Lovász: Multifractal network generator, PNAS, 2010

Installation
=============

Download the source or clone from the github repository. The letter::

    git clone git@github.com:horvatha/mfng.git

Just use the standard method::

    sudo python setup.py install

On Linux, you need to compile the file iterate.cpp if you want to use it
and put the directory of the compiled file into the PATH. I do not know
how to do it on Windows.::

    cd scripts
    make
    echo PATH=${PATH}:<path_to_the_scripts_directory> >> ~/.bashrc

Tutorial and library reference
==================================
The documentation have been included in `the documentation of the cxnet
package <https://pythonhosted.org/cxnet/mfng_en.html>`_.
