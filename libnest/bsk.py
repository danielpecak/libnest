#!/usr/bin/env python3
# -*- coding:utf-8 -*-
# =========== Description
# Original paper with formulas
# https://journals.aps.org/prc/pdf/10.1103/PhysRevC.80.065804
# NOTE: Table I is outdated; use the following data:
# forBSk31
# According to Goriely, Chamel, Pearson PRC 93 034337 (2016)
#================================
#define T0  -2302.01  // [MeV*fm^3]
#define T1  762.99    // [MeV*fm<sup>5</sup>]
#define T2  0.0       // [MeV*fm<sup>5</sup>]
#define T3  13797.83  // [MeV*fm^(3+3*ALPHA)]
#define T4  -500      // [MeV*fm^(5+3*BETA)]
#define T5  -40       // [MeV*fm^(5+3*GAMMA)]
#define X0  0.676655  // [1]
#define X1  2.658109  // [1]
#define T2X2 -422.29  // [MeV*fm<sup>5</sup>]
#define X3  0.83982   // [1]
#define X4  5.        // [1]
#define X5  -12.      // [1]
#define W0  62.174   // [MeV*fm<sup>5</sup>]
#define ALPHA (1./5.)  // [1]
#define BETA  (1./12.) // [1]
#define GAMMA (1./4.)  // [1]
#define YW  2.       // [1]
#define FNP 1.00     // [1]
#define FNM 1.06     // [1]
#define FPP 1.00     // [1]
#define FPM 1.04     // [1]
#define KAPPAN -36630.4 // [MeV*fm<sup>8</sup>]
#define KAPPAP -45207.2 // [MeV*fm<sup>8</sup>]
#================================
# TODO list:
# Formulas from https://journals.aps.org/prc/pdf/10.1103/PhysRevC.80.065804:
# Code formulas: A13, and dependent (A14, A15)
#
# TODO list:
# https://journals.aps.org/prc/pdf/10.1103/PhysRevC.82.035804
# Code function: Eq. 10 and plot it as a function of rho (I suppose this is what is done in
# Fig 5 in the inset). q=n or p (neutron or proton counterpart)
# rho = rho_n + rho+p
# Make plots for neutron matter (NeuM) where rho = rho_n [rho_p =0 ]
# Make plots for symmetric matter (SM) where rho = 2*rho_n [rho_p =rho_n ]
#
# TODO list:
# PHYSICAL REVIEW C 104, 055801 (2021)
# NOTE: densities such as TAU, MU, J will be given in the future from the data
# Formulas 9-14,23, A8-A10
