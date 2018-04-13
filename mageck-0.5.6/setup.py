#!/usr/bin/env python
'''
MAGeCK set up script
'''


from __future__ import print_function;

import os
import sys
from distutils.core import setup, Extension
from subprocess import call as subpcall
from distutils.command.install import install as DistutilsInstall

try:
   from distutils.command.build_py import build_py_2to3 as build_py
except ImportError:
   from distutils.command.build_py import build_py

exec(open('mageck/version.py').read())

def compile_rra():
  #
  os.chdir('rra');
  subpcall('make',shell=True);
  rev=subpcall('../bin/RRA',shell=True);
  os.chdir('../');
  return rev;

def compile_gsea():
  #
  os.chdir('gsea');
  subpcall('make',shell=True);
  rev=subpcall('../bin/mageckGSEA',shell=True);
  os.chdir('../');
  return rev;



class RRAInstall(DistutilsInstall):
  def run(self):
    # compile RRA
    if(compile_rra()!=0):
      print("CRITICAL: error compiling the RRA source code. Please check your c compilation environment.",file=sys.stderr);
      sys.exit(1);
    if(compile_gsea()!=0):
      print("CRITICAL: error compiling the mageckGSEA source code. Please check your c compilation environment.",file=sys.stderr);
      sys.exit(1);
    DistutilsInstall.run(self)



def main():
  # check python version
  if float(sys.version[:3])<2.7:
    sys.stderr.write("CRITICAL: Python version must be >=2.7!\n")
    sys.exit(1);

  setup(name='mageck',
    version=__version__,
    description='Model-based Analysis of Genome-wide CRISPR-Cas9 Knockout',
    author='Wei Li, Han Xu, Johannes Koster',
    author_email='li.david.wei@gmail.com',
    url='http://mageck.sourceforge.net',
    packages=['mageck'],
    scripts=['bin/mageck'],
    package_dir={'mageck':'mageck'},
    cmdclass={'install':RRAInstall, 'build_py': build_py},
    package_data={'mageck':['*.Rnw','*.RTemplate','*.gmt']},
    data_files=[('bin', ['bin/RRA','bin/mageckGSEA'])]
    #package_data={'mageck':['mageck/Makefile','mageck/src/*.c','include/*','utils/*']}
    #data_files=[('',['Makefile','src/*.c','include/*','utils/*'])]
  );


if __name__ == '__main__':
  main();
