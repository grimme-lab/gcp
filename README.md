# Geometrical Counter-Poise Correction

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
meson install -C _build
```

This might require administrator access.
You can alter the install prefix by ``meson configure _build --prefix=/path/to/install``.


## Contributors

- J. Gerit Brandenburg ([**@gbrandenburg**](https://github.com/gbrandenburg))
- Sebastian Ehlert ([**@awvwgk**](https://github.com/awvwgk))
- Holger Kruse ([**@hokru**](https://github.com/hokru))


## References

Please cite the GCP reference publication for work done with this program

1. H. Kruse, S. Grimme *J. Chem. Phys.* 136, 154101 (2012).
   DOI: [10.1063/1.3700154](https://doi.org/10.1063/1.3700154)
2. For periodic GCP also cite:
   J. G. Brandenburg, M. Alessio, B. Civalleri, M. F. Peintinger,
   T. Bredow, S.Grimme J. Phys. Chem. A 117, 9282-9292 (2013).
   DOI: [10.1021/jp406658y](https://doi.org/10.1021/jp406658y)


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
