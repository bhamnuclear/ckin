# ckin.c

Kinematics tool including penetrability calculation via inclusion of WCLBES
(and its dependencies from CERNLIBS).

This code can be used to calculate reaction and decay Q-values, print out
atomic masses and binding energies (B.E.s) and calculate the corresponding
Semi-Empirical Mass Formula (SEMF) value of the binding energy.

Two-body kinematics calculation (both relativistic and non-relativistic) is
included. Of additional interest, particularly for spectrometer measurements
and resonant reactions is the ability to solve for excitation energy, given a
particular ejectile energy.

For the Coulomb penetrability calculations, these can either be for an
individual orbital angular momentum value, *l*, and excitation energy, or for a
range of one or both, with *l* values from 0 &rarr; 8.

Allowed special notation for particles/decay modes are:
 
- *g*   &nbsp;&nbsp;gamma,
- *e-*  &nbsp;electron (same as *b-*),
- *b-*  &nbsp;beta- (same as *e-*),
- *e+*  &nbsp;positron (same as *b+*),
- *b+*  &nbsp;beta+ (same as *e+*),
- *ec*  &nbsp;electron capture,
- *ee*  &nbsp;double beta- decay (same as *bb*),
- *bb*  &nbsp;double beta- decay (same as *ee*).

For isotopes, allowed special notation is:

- *n*  &nbsp;<sup>1</sup>n - neutron,
- *p*  &nbsp;<sup>1</sup>H - hydrogen (proton with bound electron),
- *d*  &nbsp;<sup>2</sup>H - deuterium (deuteron with bound electron),
- *t*  &nbsp;<sup>3</sup>H - tritium (triton with bound electron),
- *a*  &nbsp;<sup>4</sup>He - helium-4 (alpha particle with bound electrons).

The input to ckin is the mass table from the latest (Open Access)
mass evaulation:  
*The Ame2020 atomic mass evaluation (I)*,  
W.J. Huang, M. Wang, F.G. Kondev, G. Audi and S. Naimi,  
Chinese Physics C**45**, 030002, March 2021,  
[DOI: 10.1088/1674-1137/abddb0](https://doi.org/10.1088/1674-1137/abddb0)  
and  
*The Ame2020 atomic mass evaluation (I)*,  
M. Wang, W.J. Huang, F.G. Kondev, G. Audi and S. Naimi,
Chinese Physics C**45**, 030003, March 2021,  
[DOI: 10.1088/1674-1137/abddaf](https://doi.org/10.1088/1674-1137/abddaf)

The table itself can be found at
https://www-nds.iaea.org/amdc/ame2020/mass_1.mas20.txt .

The previous evaluation table from 2016 can be found at  
https://www-nds.iaea.org/amdc/ame2016/mass16.txt .

The full list of ckin options is:

1. Enter/change a reaction or decay;
2. Calculate the reaction Q value;
3. Work out some kinematics;
4. Input ejectile energy and solve for excitation energy;
5. Mass look-up and SEMF/B.E.s;
6. Calculate resonant-reaction beam energy for a given excitation energy;
7. Calculate penetrabilities and reduced widths;
0. Quit.
