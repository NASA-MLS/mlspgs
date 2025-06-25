# MLSPGS
MLS Product Generation Software for Level 1 and Level 2.

This code was developed for the Aura [Microwave Limb Sounder] (https://mls.jpl.nasa.gov) project to produce daily radiance (Level 1) and geophysical (Level 2) data products.  

It is largely in Fortran and was originally under version control in a CVS repository, but in order to retain it's change history, we used [cvs2git](https://www.mcs.anl.gov/~jacob/cvs2svn/cvs2git.html) to handle the bulk of the history import, though some manual patching to the history was necessary along the way (the source code was untouched in this).

Please make sure to review the following documentation before trying to compile the software:

- [DEPENDENCIES.md](DEPENDENCIES.md) - information on what you'll need to get your development set up
- [MakeFC.md](MakeFC.md) - how to create a configurable Makefile
