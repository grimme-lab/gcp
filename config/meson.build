# This file is part of mctc-gcp.
# SPDX-Identifier: LGPL-3.0-or-later
#
# mctc-gcp is free software: you can redistribute it and/or modify it under
# the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# mctc-gcp is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# Lesser GNU General Public License for more details.
#
# You should have received a copy of the Lesser GNU General Public License
# along with mctc-gcp.  If not, see <https://www.gnu.org/licenses/>.

os = host_machine.system()

if os == 'windows'
  add_project_link_arguments(
    '-Wl,--allow-multiple-definition',
    language: 'fortran',
  )
endif

fc = meson.get_compiler('fortran')
fc_id = fc.get_id()

if fc_id == 'gcc'
  add_project_arguments(
    '-ffree-line-length-none',
    '-fbacktrace',
    language: 'fortran',
  )
elif fc_id == 'intel'
  add_project_arguments(
    '-traceback',
    language: 'fortran',
  )
elif fc_id == 'pgi'
  add_project_arguments(
    '-Mbackslash',
    '-Mallocatable=03',
    '-traceback',
    language: 'fortran',
  )
endif

# Create the tool chain library as subproject
mctc_prj = subproject(
  'mctc-lib',
  version: '>=0.1',
  default_options: [
    'default_library=static',
  ],
)
mctc_dep = mctc_prj.get_variable('mctc_dep')
lib_deps += mctc_dep
inc_dirs += mctc_prj.get_variable('mctc_inc')

if install
  install_data(
    mctc_prj.get_variable('mctc_lic'),
    install_dir: get_option('datadir')/'licenses'/meson.project_name()/'mctc-lib'
  )
endif
