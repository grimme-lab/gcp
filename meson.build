# This file is part of mctc-gcp.
#
# mctc-gcp is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# mctc-gcp is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with mctc-gcp.  If not, see <https://www.gnu.org/licenses/>.

project(
  'mctc-gcp',
  'fortran',
  version: '2.1.0',
  license: 'LGPL-3.0-or-later',
  meson_version: '>=0.53',
  default_options: [
    'buildtype=debugoptimized',
    'default_library=both',
  ],
)
install = not (meson.is_subproject() and get_option('default_library') == 'static')

# General configuration information
inc_dirs = []
lib_deps = []
subdir('config')

# Collect source of the project
srcs = []
subdir('src')

# Library target
gcp_lib = library(
  meson.project_name(),
  sources: srcs,
  version: meson.project_version(),
  dependencies: lib_deps,
  include_directories: inc_dirs,
  install: install,
)

# Export dependency for other projects and test suite
gcp_inc = gcp_lib.private_dir_include()
gcp_dep = declare_dependency(
  link_with: gcp_lib,
  include_directories: gcp_inc,
  dependencies: lib_deps,
)

# Add applications
subdir('app')

# Package the license files
gcp_lic = files(
  'COPYING',
)

if install
  # Distribute the license files in share/licenses/<name>
  install_data(
    gcp_lic,
    install_dir: get_option('datadir')/'licenses'/meson.project_name()
  )

  pkg = import('pkgconfig')
  pkg.generate(
    gcp_lib,
    description: 'Geometrical counterpoise correction',
  )
endif

# add the testsuite
subdir('test')
