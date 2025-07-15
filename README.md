# MLSPGS

Aura MLS Product Generation Software (MLSPGS) for Level 1 and Level 2.

This code and limited associated documentation herein (more documentation is elsewhere) has been developed for the Aura [Microwave Limb Sounder](https://mls.jpl.nasa.gov) project to produce daily calibrated radiance and engineering (Level 1) and geophysical data along the measurement footprint (Level 2) data products.

## More information on Aura MLS

Additional information can be found at:

> Waters, J.W. et al., "The Earth Observing System Microwave Limb Sounder (EOS MLS) on the Aura satellite," *IEEE Trans. Geosci. Remote Sensing* **44**, no. 5, [doi:10.1109/TGRS.2006.873771](https://dx.doi.org/doi:10.1109/TGRS.2006.873771), May 2006.

## Function and purpose of the software

The intent of this software is to, on a routine basis, process the raw telemetry ("Level 0" data) from the MLS instrument and generate higher level data products for use in atmospheric science research.  The vast majority of the MLS-based atmospheric science research derives from the "Level 2" data, which are vertical profiles of atmospheric composition and other data products along the orbit/measurement track.  These are derived from "Level 1" data, being calibrated microwave limb radiances (and supporting engineering telemetry).  Additional information on data processing levels can be found on the [NASA Earthdata website](https://www.earthdata.nasa.gov/learn/earth-observation-data-basics/data-processing-levels).

In addition to its use for routine data production, this software can be configured to perform other tasks, including those needed for testing (e.g., generating simulated Level 1 datasets from a known Level 2 dataset, pre-computing some look-up-table-like terms used in the Level 2 processing).  In addition to the Level 1 and Level 2 processing software, smaller programs provide additional support.

## Implementation and operational usage of the software

The software has been run on a routine basis over the 20-plus-years of MLS operations on the MLS Science-Investigator-led Processing System (SIPS).  The Level 2 processing is a notably intensive computation, requiring of order eight hours to process 24-hours of data on a ~350-cpu cluster.

The vast majority of the software is written in Fortran, requiring a Fortran-2000 compliant compiler.  It is known to work on Intel-based machines using the Intel Fortran compiler, though, in principal, should be readily compiled and run on other architectures.

Level 1 processing is a relatively lightweight task that can be run in a short order for 24-hours of data on a single CPU machine.  The Level 1 data are the largest datasets for the mission (though not particularly large by today's standards) requiring around 4Gb per day of data.

Level 2 is run in an "embarassingly parallel" manner, with the (admittedly dated) PVM (Parallel Virtual Machine) message passing system used for the (relatively lightweight) communication between the tasks.  Fundamentally, each CPU runs a single-threaded task to process a fraction (~15 degrees along the orbit track, about 4-5 minutes-worth of data).

## History of the software and this repository

The software is largely in Fortran and was originally under version control in a CVS repository, but in order to retain it's change history, we used [cvs2git](https://www.mcs.anl.gov/~jacob/cvs2svn/cvs2git.html) to handle the bulk of the history import, though some manual patching to the history for a large number of files was necessary along the way (the source code was untouched in this).  The original CVS repository, consisting of a heirarchy of RCS ("...,v") files is stored untouched in the `original-cvs-repository-23jun2025` branch.  This branch also includes a copy (taken out of context) of the (homespun) tools and scripts to patch the affected files and perform the cvs2-git importing.

### Branches

- `main`:  branch with the official release of the software
- `njl-cvs-import-20jun25`: branch as imported by `cvs2git`
- `original-cvs-repository-23jun2025`: original CVS repository

## Other files within this repository

Please make sure to review the following documentation before trying to compile the software:

- [DEPENDENCIES.md](DEPENDENCIES.md) - information on what you'll need to get your development set up
- [MakeFC.md](MakeFC.md) - how to create a configurable Makefile
