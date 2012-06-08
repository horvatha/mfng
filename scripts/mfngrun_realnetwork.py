#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
from __future__ import print_function
import cmfng
import cxnet
from cxnet.archives import get_netdata_directory
import os

filename = "as-"
results = cxnet.archives.get_archive_name(filename)
archive, file_ = results[0]
file_name = os.path.join(get_netdata_directory(), "netdata", "{0}.gml".format(file_))
print("file name = {0}, archive = {1}".format(file_name, archive))
degrees = cxnet.Graph.Read_GML(file_name).degree()
n = len(degrees)
print("archive: {0}, n = {1}, max_degree = {2}".format(archive, n, max(degrees)))

for steps in [1000]*2:
    T0 = 0.2
    Tlimit = T0/10000
    generator = cmfng.Generator(T0=T0, steps=steps, Tlimit=Tlimit,
            m=3, K=4,
            n=n,
            divexponent = 7,
            project = "2012_jan_realnetworks",
            )
    generator.append_property(
            cmfng.LogBinnedDegDist(
                degrees = degrees,
                )
            )
    generator.go()

