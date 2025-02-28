#! /usr/bin/env python

# piPipes, a set of pipelines for PIWI-interacting RNA (piRNA) and transposon analysis
# Copyright (C) 2014  Bo Han, Wei Wang, Zhiping Weng, Phillip Zamore
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

import sys
from collections import defaultdict #, Counter # Counter is not available in python until 2.7
ct = defaultdict ( lambda: defaultdict (int) )
ext = int (sys.argv[1]) * 2 + 1
for line in sys.stdin:
        seq = str (line.split()[1])
        for i in range (0,ext):
                ct[i][seq[i].upper()]+=1
for n in range (0, ext):
        A = ct[n]['A']
        C = ct[n]['C']
        G = ct[n]['G']
        T = ct[n]['T']
        print(str (A) + '\t' + str (C) + '\t' + str (G) + '\t' + str (T))
