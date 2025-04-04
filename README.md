# ckin
Kinematics tool including penetrability calculation via inclusion of WCLBES (and its dependencies from CERNLIBS).

This code can be used to calculate reaction and decay Q-values, print out atomic masses and binding energies (B.E.s) and calculate the corresponding Semi-Empirical Mass Formula (SEMF) value of the binding energy.

Two-body kinematics calculation (both relativistic and non-relativistic) is included. Of additional interest, particularly for spectrometer measurements and resonant reactions is the ability to solve for excitation energy, given a particular ejectile energy.

For the Coulomb penetrability calculations, these can either be for an individual orbital angular momentum value, l, and excitation energy, or for a range of one or both, with l values from 0 --> 8.

Allowed special notation for particles/decay modes are:  g, e-, b-, e+, b+, ec, ee, bb.
For isotopes, allowed special notation is:  n (1n), p (1H), d (2H), t (3H), a (4He).

The input to ckin is the mass table from the latest mass evaulation
(e.g. AME 2020 by W.J. Huang, M. Wang, F.G. Kondev, G. Audi and S. Naimi, Chinese Physics C45, 030002, March 2021).

The full list of ckin options is:
 1) Enter/change a reaction or decay.
 2) Calculate the reaction Q value.
 3) Work out some kinematics.
 4) Input ejectile energy and solve for excitation energy.
 5) Mass look-up and SEMF/B.E.s.
 6) Calculate resonant-reaction beam energy for a given excitation energy.
 7) Calculate penetrabilities and reduced widths.

