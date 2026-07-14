# SPDX-License-Identifier: MIT
# Copyright (c) 2022-2026 Daniel Pęcak
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
        Works for 1D, 2D, and 3D grids. A dimension not listed in the file
        defaults to 1 (a flat direction).

    Args:
        PATH (string): location path
        PREFIX (string): prefix of the file

    Returns:
        list of ints: a list of points in each direction [NX, NY, NZ].
    """
    file = PATH + PREFIX + '_info.txt'
    with open(file) as f:
        lines = f.readlines()
    NX = NY = NZ = 1  # missing dimensions default to 1 (flat direction)
    for l in lines:
        if 'NX' in l: NX = int(l.split('=')[-1])
        if 'NY' in l: NY = int(l.split('=')[-1])
        if 'NZ' in l: NZ = int(l.split('=')[-1])
    return [NX, NY, NZ]

def txt2df(file,sufix,cols):
    r"""Reads \*.txt files to pandas dataframe.

    Args:
        file (string): path and prefix to the data
        sufix (string): sufix pointing to the type of data e. g. 'density', 'delta' etc.
        cols (list of strings): names of columns in the file

    Returns:
        pandas.DataFrame: data about spatial extend of observables

    """
    file = file + '_' + sufix + '.txt'
    df = pd.read_csv(file, comment='#', sep=r'\s+', header=None)
    df = df.set_axis(cols, axis=1)  # pandas >= 2.0: set_axis is not in-place
    return df



if __name__ == '__main__':
    pass
