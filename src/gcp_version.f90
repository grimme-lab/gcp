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

!> Versioning information on this library.
module gcp_version
   implicit none
   private

   public :: gcp_version_string, gcp_version_compact
   public :: get_gcp_version


   !> String representation of the mctc-gcp version
   character(len=*), parameter :: gcp_version_string = "2.3.1"

   !> Numeric representation of the mctc-gcp version
   integer, parameter :: gcp_version_compact(3) = [2, 3, 1]


contains


!> Getter function to retrieve mctc-gcp version
subroutine get_gcp_version(major, minor, patch, string)

   !> Major version number of the mctc-gcp version
   integer, intent(out), optional :: major

   !> Minor version number of the mctc-gcp version
   integer, intent(out), optional :: minor

   !> Patch version number of the mctc-gcp version
   integer, intent(out), optional :: patch

   !> String representation of the mctc-gcp version
   character(len=:), allocatable, intent(out), optional :: string

   if (present(major)) then
      major = gcp_version_compact(1)
   end if
   if (present(minor)) then
      minor = gcp_version_compact(2)
   end if
   if (present(patch)) then
      patch = gcp_version_compact(3)
   end if
   if (present(string)) then
      string = gcp_version_string
   end if

end subroutine get_gcp_version


end module gcp_version
