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
# $ ./units.py
import numpy as np

DENSEPSILON = 1e-12

HBARC=197.         # \hbar c [MeV fm]
VUNIT=197./940*100 # alpha constant over the neutron mass (100 to have percentage)
hbar22M0   =20.72  # neutron bare mass
MN   =939.56542052 # neutron mass [MeV/c]
MP   =938.27208816 # proton mass [MeV/c]
T0   =-2302.01     # skyrme parameter t0 [MeV*fm^3]
T1   =762.99       # skyrme parameter t1 [MeV*fm<sup>5</sup>]
T2   =0.0          # skyrme parameter t2 [MeV*fm<sup>5</sup>]
T3   =13797.83     # skyrme parameter t3 [MeV*fm^(3+3*ALPHA)]
T4   =-500         # skyrme parameter t4 [MeV*fm^(5+3*BETA)]
T5   =-40          # skyrme parameter t5 [MeV*fm^(5+3*GAMMA)]
X0   =0.676655     # skyrme parameter x0 [1]
X1   =2.658109     # skyrme parameter x1 [1]
T2X2 =-422.29      # skyrme parameter x2t2 [1][MeV*fm<sup>5</sup>]
X3   =0.83982      # skyrme parameter x3 [1]
X4   =5.           # skyrme parameter x4 [1]
X5   =-12.         # skyrme parameter x5 [1]
ALPHA =(1./5.)     # [1]
BETA =(1./12.)     # [1]
GAMMA =(1./4.)     # [1]
YW   =2.           # [1]
FNP  =1.00         # [1]
FNM  =1.06         # [1]
FPP  =1.00         # [1]
FPM  =1.04         # [1]
KAPPAN =-36630.4   # [MeV*fm<sup>8</sup>]
KAPPAP =-45207.2   # [MeV*fm<sup>8</sup>]

def printhbar22M0():
    """Prints hbar22M0 constant"""
    print(hbar22M0)

def KtoMev(temp):
    """Function converts temperature units from Kelvins to MeVs"""
    return temp/11604525006.1598

def MeVtoK(temp):
    """Function converts temperature units from MeVs to Kelvins"""
    return temp*11604525006.1598

def fm3togcm3(rho):
    """Function converts units 1/fm^3 into g/cm^3."""
    return rho*1.67e15

def gcm3tofm3(rho):
    """Function converts units g/cm^3 into 1/fm^3."""
    return rho/1.67e15

if __name__ == '__main__':
    print("# Unit conversion: 1/fm^3 into g/cm^3")
    print("# 1/fm^3 = 1.67 * 10^15 g/cm^3")
    print("{0:.2E}".format(fm3togcm3(1.)))
    print("{0:.2E}".format(fm3togcm3(.08)))
    print("{0:.2E}".format(fm3togcm3(.05)))
    print("{0:.2E}".format(fm3togcm3(.03)))

    print("\n\n")
    print("# Unit conversion: g/cm^3 into 1/fm^3.")
    print("# g/cm^3 = 0.59 * 10^-15 1/fm^3")
    print("# Nuclear saturation: {0:.2E}".format(gcm3tofm3(3.e14)))
    print("# Neutron drip: {0:.2E}".format(gcm3tofm3(4.e11)))
    print("{:.6f}".format(gcm3tofm3(4.e11)))

rhoSAT = gcm3tofm3(3.e14)
rhoND  = gcm3tofm3(4.e11)
