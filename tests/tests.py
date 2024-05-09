#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# =========== Info: who, where, when
# @author Daniel Pęcak <Daniel.Pecak@pw.edu.pl>
# Warsaw Technical University, Université Libre de Bruxelles
# On leave: Institute of Physics, Polish Academy of Sciences, Warsaw
# March 2022, Brussels
# =========== Description
# Script for converting the density units for nuclear matter
# in the ranges typical for neutron stars.
# =========== Usage example
# $ ./units_test.py


from libnest.units import *
import numpy as np

mylist = [0,1,2,3,4]
print(mylist)
nplist = np.array(mylist)
print(nplist)

print(KtoMev(mylist))
print(KtoMev(nplist))
print(MeVtoK(mylist))
print(MeVtoK(nplist))
print(fm3togcm3(mylist))
print(fm3togcm3(nplist))
print(gcm3tofm3(mylist))
print(gcm3tofm3(nplist))

if __name__ == '__main__':
    pass
