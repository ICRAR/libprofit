#
#    ICRAR - International Centre for Radio Astronomy Research
#    (c) UWA - The University of Western Australia, 2016
#    Copyright by UWA (in the framework of the ICRAR)
#    All rights reserved
#
#    This library is free software; you can redistribute it and/or
#    modify it under the terms of the GNU Lesser General Public
#    License as published by the Free Software Foundation; either
#    version 2.1 of the License, or (at your option) any later version.
#
#    This library is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public
#    License along with this library; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place, Suite 330, Boston,
#    MA 02111-1307  USA
"""
A small tool to run several
"""

import re
import shutil
import os
import subprocess
import itertools
import collections
import argparse
import multiprocessing

compilers = {
    'g++': (4.9, 6),
    'clang++': (3.6, 3.7, 3.8),
}

settings = {
    'osize': '-g -Os -fPIC',
    'o1':    '-g -O1',
    'o2':    '-g -O2',
    'o3':    '-g -O3',
    'o2n':   '-g -O2 -march=native',
    'o3n':   '-g -O3 -march=native',
}

config = collections.namedtuple('config', 'comp version name cxxflags')
def get_configs():
    versioned_comps = [(comp,v) for comp,versions in compilers.items() for v in versions]
    all_configs = itertools.product(versioned_comps, settings.items())
    return [config(c[0][0], c[0][1], c[1][0], c[1][1]) for c in all_configs]


root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

def _cmake_c(c):
    dirname = os.path.join(root, 'builds', c.comp, str(c.version), c.name)
    if not os.path.isdir(dirname):
        os.makedirs(dirname)

    cmake_cxxflags = '-DCMAKE_CXX_FLAGS=%s' % c.cxxflags
    cmake_cxx = '-DCMAKE_CXX_COMPILER=%s-%s' % (c.comp, str(c.version))
    subprocess.check_call(['rm', '-rf'], cwd=dirname)
    subprocess.check_call(['cmake', cmake_cxx, cmake_cxxflags, root], cwd=dirname)

def cmake(jobs):
    bdir = os.path.join(root, 'builds')
    if not os.path.isdir(bdir):
        os.mkdir(bdir)
    p = multiprocessing.Pool(jobs)
    p.map(_cmake_c, get_configs())
    p.close()

def _make_c(c):
    dirname = os.path.join(root, 'builds', c.comp, str(c.version), c.name)
    subprocess.check_call(['make'], cwd=dirname)

def make(jobs):
    p = multiprocessing.Pool(jobs)
    p.map(_make_c, get_configs())
    p.close()

_r = re.compile('Ran [0-9]+ iterations in [0-9\.]+ \[s\] \(([0-9\.]+) \[ms\] per iteration\)')
def benchmark(args=None):
    for c in get_configs():
        dirname = os.path.join('builds', c.comp, str(c.version), c.name)
        cmdline = ['./profit-cli']
        if args:
            cmdline += args
        else:
            cmdline += ['-w', '400', '-H', '400', ]
            cmdline += ['-p', 'sersic:xcen=201:ycen=201:mag=15.87:re=4.65:nser=4.7:ang=-21:axrat=0.74:box=0']
            cmdline += ['-p', 'sersic:xcen=201:ycen=201:mag=16.63:re=41:nser=1:axrat=0.29:box=0']
            cmdline += ['-i', '100']
        out = subprocess.check_output(cmdline, cwd=dirname).strip()
        ms_per_iter = _r.match(out).group(1)
        print('%s, %s, %s, %s' % (c.comp, str(c.version), c.name, ms_per_iter))

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="A tool to run commands through different builds")
    g = parser.add_mutually_exclusive_group(required=True)
    g.add_argument('-c', '--cmake', help="Runs CMake to prepare the build", action="store_true")
    g.add_argument('-m', '--make', help="Runs make in the builds", action="store_true")
    g.add_argument('-b', '--benchmark', help="Runs a benchmark in the builds. -j has no effect", action="store_true")
    parser.add_argument('-j', help="Numbers of jobs to run in parallel (processes)", type=int, dest="jobs", default=1)
    parser.add_argument('args', help="Arguments to use for benchmarking (optional)", nargs=argparse.REMAINDER)

    args = parser.parse_args()
    if args.cmake:
        cmake(args.jobs)
    if args.make:
        make(args.jobs)
    elif args.benchmark:
        # -- is there? it shouldn't...
        if len(args.args) > 0 and args.args[0] == '--':
            args.args = args.args[1:]
        benchmark(args.args)