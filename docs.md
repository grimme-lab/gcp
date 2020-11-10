---
project: mctc-gcp
summary: Geometrical Counterpoise Correction
project_github: https://github.com/grimme-lab/gcp
project_download: https://github.com/grimme-lab/gcp/releases
author: Grimme group, Bonn
github: https://github.com/grimme-lab
src_dir: ./src
output_dir: ./docs
exclude_dir: ./test
docmark: <
predocmark: >
source: true
graph: false
sort: alpha
print_creation_date: true
extra_mods: iso_fortran_env:https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fFORTRAN_005fENV.html
creation_date: %Y-%m-%d %H:%M %z
md_extensions: markdown.extensions.toc
               markdown.extensions.smarty
---

[![DOI](https://img.shields.io/badge/DOI-10.1063%2F1.3700154-blue)](https://doi.org/10.1063/1.3700154)
[![DOI](https://img.shields.io/badge/DOI-10.1021%2Fjp406658y-blue)](https://doi.org/10.1021/jp406658y)

Program to compute geometrical counterpoise corrections.


## Getting Started

### Meson

Create a new meson project and include `gcp` either as git-submodule in your subprojects directory or create a wrap file to fetch it from upstream:

```ini
[wrap-git]
directory = mctc-gcp
url = https://github.com/grimme-lab/gcp
revision = head
```

The `gcp` library depends on the [MCTC-library](https://github.com/grimme-lab/mctc-lib).
You might have to promote those dependencies to your subprojects by copying the wrap files from the `gcp` subprojects.

To load the project the necessary boilerplate code for subprojects is just

<!--pygments doesn't know about meson, python highlighting looks okayish-->
```python
gcp_prj = subproject(
  'mctc-gcp',
  version: '>=0.1',
  default_options: [
    'default_library=static',
  ],
)
gcp_dep = gcp_prj.get_variable('gcp_dep')
```

Now you can add `gcp_dep` to your dependencies and access the public API by the `gcp` module.

We recommend to set the default library type of `gcp` to static when linking your applications or library against it.
Note for library type both and shared `gcp` will install itself along with your project.

For more fine-tuned control you can access:

- the library target with `gcp_lib`
- the private include dir of this target, containing the Fortran module files, with `gcp_inc`
- the license files of `gcp` with `gcp_lic`

If you are linking your application statically against `gcp` and still want to distribute the license files of `gcp` (thank you), just use

```python
install_data(
  gcp_prj.get_variable('gcp_lic'),
  install_dir: get_option('datadir')/'licenses'/meson.project_name()/'mctc-gcp',
)
```


### Fortran Package Manager (fpm)

This project supports [fpm](https://github.com/fortran-lang/fpm) as build system as well.
Just add it to the dependencies in your `fpm.toml` file:

```toml
[dependencies]
[dependencies.mctc-gcp]
git = "https://github.com/grimme-lab/gcp"
```
