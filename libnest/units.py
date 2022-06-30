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



HBARC=197.         # \hbar c [MeV fm]
VUNIT=197./940*100 # alpha constant over the neutron mass (100 to have percentage)
M0   =20.72        # neutron bare mass

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
