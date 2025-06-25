# Installation Dependencies

## NOTE

It is possible that you won't need to read this document. If you obtained the
mls software as part of a shrinkwrap package, look there for the files
build-xxx (where xxx might be "nag"), env.linux, and f90. After judiciously
editing these, just cd to their directory and type `build-xxx`.

## SUMMARY

Before building the mls software you will need the following

 1. f95 compiler (NAG, Lahey, or Sun)
 2. Installed SDP Toolkit libraries:
    1. PGSTK
    2. HDFEOS
    3. HDF
 3. mls software (probably where you're reading this)
 4. PVM--a library to emulate parallel virtual machines
 5. FFTW--a library to accomplish fast Fourier tranforms
 6. BLAS and LAPACK--libraries to speed up linear algebra

The first step is to find out if you already have (1)-(6). If you do, you may
safely ignore the rest of this guide. Obtaining (1), (3), and (6) is not
within the scope of this guide, although if you are reading it online, you
probably already have (3).

The rest of this guide focuses on how to obtain (2), (5), and (6). After reading
it you will probably want to read MakeFC.md or else type
> make -f MakeFC disthelp

## THE SDP TOOLKIT

The distribution to build the libraries and ancillary files making up the SDP
Toolkit can be found at:

> [https://edhs1.gsfc.nasa.gov/](https://edhs1.gsfc.nasa.gov/)

It comes with instructions for building it for a number of supported platforms.
Additional information on how to build it for the linux OS can be found at

> [https://www.met.ed.ac.uk/~hcp/hdfeos_notes.txt](https://www.met.ed.ac.uk/~hcp/hdfeos_notes.txt)

After building this successfully, keep careful note of its location. An easy
way to do this is to source pgs-env.csh or pgs-dev-env.csh in your login shell,
assuming you use the csh or tcsh shells. That way automatically defines three
environment variables that will make configuring the mls software easier:

> HDFEOS_LIB
>
> HDFLIB
>
> PGSLIB

If you get the latest version of the toolkit, you will automatically build two
added libraries: hdf5 and hdfeos5. These too will be needed for correctly
building the mls software if you use the latest version.

## PVM

The distribution to build the libraries and ancillary files making up the PVM
software system can be found at

> [https://www.netlib.org/pvm3/index.html](https://www.netlib.org/pvm3/index.html)

in particular the file: [https://www.netlib.org/pvm3/pvm3.4.3.tgz](https://www.netlib.org/pvm3/pvm3.4.3.tgz)

It comes with instructions for building for a unix-style platforms. In brief,
after downloading the file, say to your home directory, issue the following
commands

```sh
gunzip pvm3.4.3.tgz
tar -xvf pvm3.4.3.tar
cd pvm3
setenv PVM_ROOT `pwd`
make
```

Note carefully the value of `PVM_ROOT` --it will be needed when configuring the
mls software.

## FFTW

It is possible that workable fftw libraries already exist on your platform. For
example, JPL's EOS-MLS computers all have it in `/usr/local/lib`. You can check
this via

`me:51% ls /usr/local/lib/\*fftw*`

> /usr/local/lib/libfftw.a
>
> /usr/local/lib/librfftw.a
>
> /usr/local/lib/libfftw.la*
>
> /usr/local/lib/librfftw.la*

and the result, displayed beneath the ls command, indicates a positive match. In
such a case, note that the `FFTW_ROOT` should take the value of `/usr/local/lib`.

Otherwise, the distribution to build the libraries and ancillary files making up
the FFTW libraries can be found at

> [https://www.fftw.org/download.html](https://www.fftw.org/download.html)

Follow the instructions found in `INSTALL` and note carefully where it places
libfftw.a and librfftw.a. If you lack root privileges, remember to use

```sh
configure --prefix=installation_directory_name
```

before

```sh
make
make install
```

An extra wrinkle, especially if you are installing onto a linux platform, is
that you may have to edit the `configure.in` file. Without this change, the
resulting library will contain object files with the wrong number of
underscores appended. In particular, look for the line with

> AC_CHECK_PROGS(F77, f77 xlf xlf77 cf77 fl32 g77 fort77 f90 xlf90)

and replace it with

> AC_CHECK_PROGS(F77, f95 f77 xlf xlf77 cf77 fl32 g77 fort77 f90 xlf90)

Thanks to Robert Thurstans (JPL) for pointing this out.

## BLAS

It is possible that highly-optimized blas libraries exist on your platform or
are made available through your f95 compiler or hardware vendor. If so, you
should note their paths and names. For example, Intel supplies: `libmkl32_p3.a`,
`libmkl32_p4.a` and `libmkl64_itp.a`. Or your compiler vendor or OS vendor may have
put a copy in an easy-to-access place: e.g.

`me/mlspgs:186%ls /usr/lib/\*blas*`

> /usr/lib/libblas.a
>
> /usr/lib/libblas.so.3@
>
> /usr/lib/libblas.so.3.0.3*
>
> /usr/lib/libblas.so@
>
> /usr/lib/libblas.so.3.0@

In either case, note

1. Whether you have such a library; and
2. If so, its path and name

You have to choose one of the three options:

1. Use the externally-supplied library of BLAS routines
2. Use f95 intrinsics to satisfy BLAS routines
3. Let the mls software build its own BLAS routines

Note that (2) and (3) will always be available, and (3) is the built-in default

## LAPACK

Highly-optimized LAPACK libraries may also exist on your platform. Again,
note their paths and names.

If you are ambitious, you might check the site named ATLAS for "Automatically
Tuned Linear Algebra Software"

> [https://www.netlib.org/atlas/index.html](https://www.netlib.org/atlas/index.html)

so that you can download or indeed build your own. In an albeit artificial
benchmark its use resulted in an order-of-magnitude speedup on a dual-processor
pentium III doing matrix operations.

You have to choose one of the two options:

1. Use the externally-supplied library of LAPACK routines
2. Use f95 intrinsics to satisfy LAPACK routines

Note that (2) will always be available, and is the built-in default.