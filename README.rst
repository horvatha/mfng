======================
cmfng Python package
======================

A realization of Multifractal Network Generator.

See:
Palla - Vicsek - Lovász: Multifractal network generator, PNAS, 2010

Installation
=============

Just use the standard method::

    sudo python setup.py install

On Linux, you need to compile the file iterate.cpp if you want to use
it and put the compiled file into the PATH. I do not know how to do it
on Windows.::

    cd scripts
    g++ -o iterate iterate.cpp
    echo >> ~/.bashrc PATH=${PATH}:<path_to_the_script_directory>


