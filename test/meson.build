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

tests = [
  'gcp',
]

test_srcs = files(
  'main.f90',
  'testsuite_structures.f90',
)
foreach t : tests
  test_srcs += files('test_@0@.f90'.format(t.underscorify()))
endforeach

tester = executable(
  'tester',
  sources: test_srcs,
  dependencies: gcp_dep,
)

foreach t : tests
  test(t, tester, args: t)
endforeach
