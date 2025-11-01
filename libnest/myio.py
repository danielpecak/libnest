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
"""
Module: Myio
============
This module provides functions to read and write data in different formats.

Currently supported formats:
- Text files (.txt) with columnar data
- Grid information files (_info.txt)

Planned features:
- WDATA format support using Forbes' python library

List of functions
-----------------
"""
import pandas as pd

def dimension(nx,ny,nz):
    """
    Calculates the dimensionality of the problem.

    Args:
        nx (int): number of grid points along x-axis
        ny (int): number of grid points along y-axis
        nz (int): number of grid points along z-axis

    Returns:
        int: dimensionality of the problem
    """
    flats=0
    if nx==1: flats+=1
    if ny==1: flats+=1
    if nz==1: flats+=1
    return 3-flats

def readDimTxt(PATH, PREFIX):
    """Reads data from the file ``PREFIX_info.txt`` in location ``PATH``.
    It gets information about the grid size and returns them in a list.

    .. note::
        This works only for 2D!

    Args:
        PATH (string): location path
        PREFIX (string): prefix of the file

    Returns:
        list of ints: a list of points in each direction [NX, NY, NZ].
    """
    print("# Directory:\t{} \n# Prefix:\t{}\n".format(PATH, PREFIX))
    file=PATH+PREFIX+'_info.txt'
    print("# Reading {}".format(PREFIX+'_info.txt'))
    with open(file) as f:
        lines = f.readlines()
    for l in lines:
        if 'NX' in l: NX=int(l.split('=')[-1])
        if 'NY' in l: NY=int(l.split('=')[-1])
        if 'NZ' in l: NZ=int(l.split('=')[-1])
        if 'DX' in l: DX=float(l.split('=')[-1])
        if 'DY' in l: DY=float(l.split('=')[-1])
        if 'DZ' in l: DZ=float(l.split('=')[-1])
    if (dimension(NX,NY,NZ)!=2):
        sys.exit("# The system is {}D, not 2D!".format(dimension(NX,NY,NZ)))
    if (NX!=NY) : print("# Warning: NX != NY   !")
    return([NX,NY,NZ])

def txt2df(file,sufix,cols):
    """Reads \*.txt files to pandas dataframe.

    Args:
        file (string): path and prefix to the data
        sufix (string): sufix pointing to the type of data e. g. 'density', 'delta' etc.
        cols (list of strings): names of columns in the file

    Returns:
        pandas.DataFrame: data about spatial extend of observables

    """
    file=file+'_'+sufix+'.txt'
    print("# Reading {}".format(file))
    df = pd.read_csv(file, comment='#', sep='\s+', header=None)
    df.set_axis(cols, axis=1,inplace=True)
    return df



if __name__ == '__main__':
    pass
