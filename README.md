# Geometrical Counter-Poise Correction

[![GPL-3.0-or-later](https://img.shields.io/github/license/grimme-lab/gcp)](LICENSE)
[![CI](https://github.com/grimme-lab/gcp/workflows/CI/badge.svg)](https://github.com/grimme-lab/gcp/actions)
[![DOI](https://img.shields.io/badge/DOI-10.1063%2F1.3700154-blue)](https://doi.org/10.1063/1.3700154)
[![DOI](https://img.shields.io/badge/DOI-10.1021%2Fjp406658y-blue)](https://doi.org/10.1021/jp406658y)


## Installation

To build this project from the source code in this repository you need to have
- a Fortran compiler supporting Fortran 2008
- [meson](https://mesonbuild.com) version 0.53 or newer
- a build-system backend, *i.e.* [ninja](https://ninja-build.org) version 1.7 or newer

Setup a build with

```
meson setup _build
```

You can select the Fortran compiler by the `FC` environment variable, currently this project supports GCC and Intel compilers.
To compile the project run

```
meson compile -C _build
```

You can run the projects testsuite with

```
meson test -C _build --print-errorlogs
```

If the testsuite passes you can install with

```
meson configure _build --prefix=/path/to/install
meson install -C _build
```

This might require administrator access depending on the chosen install prefix.
Now you are ready to use ``mctc-gcp``.


## Usage

To calculate the geometrical counter poise correction use the [``mctc-gcp(1)``](man/mctc-gcp.1.adoc) program.
For a Hartree-Fock calculation in split valence basis use

```
mctc-gcp coord -l hf/def2-svp
```

Similarly, other methods can be selected.
Special levels are the “3c” methods ``hf-3c``, ``pbeh-3c``, ``hse-3c``, ``b97-3c`` and ``r2scan-3c``.

Periodic calculations can be performed by providing periodic input, like Vasp's POSCAR, riper coord file or supercell DFTB+ general format.

```
mctc-gcp POSCAR -l hse-3c
```

For more details look up [the manual page](man/mctc-gcp.1.adoc).


## Contributors

- J. Gerit Brandenburg ([@gbrandenburg](https://github.com/gbrandenburg))
- Sebastian Ehlert ([@awvwgk](https://github.com/awvwgk))
- Holger Kruse ([@hokru](https://github.com/hokru))


## References

Please cite the GCP reference publication for work done with this program

1. H. Kruse, S. Grimme *J. Chem. Phys.* **136**, 154101 (2012).
   DOI: [10.1063/1.3700154](https://doi.org/10.1063/1.3700154)
2. For periodic GCP also cite:
   J. G. Brandenburg, M. Alessio, B. Civalleri, M. F. Peintinger,
   T. Bredow, S. Grimme *J. Phys. Chem. A* **117**, 9282–9292 (2013).
   DOI: [10.1021/jp406658y](https://doi.org/10.1021/jp406658y)

For the “3c” methods see:

1. R. Sure, S. Grimme *J. Comput. Chem.* **34**, 1672–1685 (2013).
   DOI: [10.1002/jcc.23317](https://doi.org/10.1002/jcc.23317)
2. S. Grimme, J. G. Brandenburg, C. Bannwarth, A. Hansen *J. Chem. Phys.* **143**,
   054107 (2015). DOI: [10.1063/1.4927476](https://doi.org/10.1063/1.4927476)
3. J. G. Brandenburg, E. Caldeweyher, S. Grimme. *Phys. Chem. Chem. Phys.* **18**,
   15519–15523 (2016).
   DOI: [10.1039/C6CP01697A](https://doi.org/10.1039/C6CP01697A)
4. J. G. Brandenburg, C. Bannwarth, A. Hansen, S. Grimme *J. Chem. Phys.* **148**,
   064104 (2018). DOI: [10.1063/1.5012601](https://doi.org/10.1063/1.5012601)
5. S. Grimme, A. Hansen, S. Ehlert, J.-M. Mewes, *ChemRxiv*, preprint (2020).
   DOI: [10.26434/chemrxiv.13333520.v1](https://doi.org/10.26434/chemrxiv.13333520.v1)


## License

This project is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This project is distributed in the hope that it will be useful,
but without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose.  See the
GNU General Public License for more details.

Unless you explicitly state otherwise, any contribution intentionally
submitted for inclusion in this project by you, as defined in the
GNU General Public license, shall be licensed as above, without any
additional terms or conditions.
