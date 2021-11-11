! This file is part of mctc-gcp.
! SPDX-Identifier: LGPL-3.0-or-later
!
! mctc-gcp is free software: you can redistribute it and/or modify it under
! the terms of the Lesser GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! mctc-gcp is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! Lesser GNU General Public License for more details.
!
! You should have received a copy of the Lesser GNU General Public License
! along with mctc-gcp.  If not, see <https://www.gnu.org/licenses/>.

!> Driver for unit testing
program tester
   use, intrinsic :: iso_fortran_env, only : error_unit
   use mctc_env_testing, only : run_testsuite, new_testsuite, testsuite_type, &
      & select_suite, run_selected
   use test_gcp_gradient, only : collect_gcp_gradient
   use test_gcp_hf_methods, only : collect_gcp_hf_methods
   use test_gcp_dft_methods, only : collect_gcp_dft_methods
   use test_gcp_3c_methods, only : collect_gcp_3c_methods
   use test_gcp_pbc, only : collect_gcp_pbc
   implicit none
   integer :: stat, is
   character(len=:), allocatable :: suite_name, test_name
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   stat = 0

   testsuites = [ &
      & new_testsuite("gcp-gradient", collect_gcp_gradient), &
      & new_testsuite("gcp-hf-methods", collect_gcp_hf_methods), &
      & new_testsuite("gcp-dft-methods", collect_gcp_dft_methods), &
      & new_testsuite("gcp-3c-methods", collect_gcp_3c_methods), &
      & new_testsuite("gcp-pbc", collect_gcp_pbc) &
      & ]

   call get_argument(1, suite_name)
   call get_argument(2, test_name)

   if (allocated(suite_name)) then
      is = select_suite(testsuites, suite_name)
      if (is > 0 .and. is <= size(testsuites)) then
         if (allocated(test_name)) then
            write(error_unit, fmt) "Suite:", testsuites(is)%name
            call run_selected(testsuites(is)%collect, test_name, error_unit, stat)
            if (stat < 0) then
               error stop 1
            end if
         else
            write(error_unit, fmt) "Testing:", testsuites(is)%name
            call run_testsuite(testsuites(is)%collect, error_unit, stat)
         end if
      else
         write(error_unit, fmt) "Available testsuites"
         do is = 1, size(testsuites)
            write(error_unit, fmt) "-", testsuites(is)%name
         end do
         error stop 1
      end if
   else
      do is = 1, size(testsuites)
         write(error_unit, fmt) "Testing:", testsuites(is)%name
         call run_testsuite(testsuites(is)%collect, error_unit, stat)
      end do
   end if

   if (stat > 0) then
      write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop 1
   end if


contains


!> Obtain the command line argument at a given index
subroutine get_argument(idx, arg)

   !> Index of command line argument, range [0:command_argument_count()]
   integer, intent(in) :: idx

   !> Command line argument
   character(len=:), allocatable, intent(out) :: arg

   integer :: length, stat

   call get_command_argument(idx, length=length, status=stat)
   if (stat /= 0) then
      return
   endif

   allocate(character(len=length) :: arg, stat=stat)
   if (stat /= 0) then
      return
   endif

   if (length > 0) then
      call get_command_argument(idx, arg, status=stat)
      if (stat /= 0) then
         deallocate(arg)
         return
      end if
   end if

end subroutine get_argument


end program tester
