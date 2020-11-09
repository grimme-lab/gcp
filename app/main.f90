program main
   use, intrinsic :: iso_fortran_env, only : output_unit, error_unit, input_unit
   use mctc_env
   use mctc_io
   use gcp, only : gcp_call, wregrad_tm, get_gcp_version
   use strings, only : lowercase
   implicit none
   character(len=*), parameter :: prog_name = "mctc-gcp"

   character(len=:), allocatable :: input
   integer, allocatable :: input_format
   type(structure_type) :: mol
   type(error_type), allocatable :: error
   real(wp) :: energy, gradlatt(3, 3), lattice(3, 3)
   real(wp), allocatable :: gradient(:, :)
   character(len=:), allocatable :: method
   character(len=20) :: lc_method
   logical :: dograd, dohess, echo, parfile
   integer, allocatable :: ifrez(:)

   call get_arguments(input, input_format, method, dograd, dohess, echo, &
      & parfile, error)
   if (allocated(error)) then
      write(error_unit, '(a)') error%message
      error stop
   end if

   if (input == "-") then
      if (.not.allocated(input_format)) input_format = filetype%xyz
      call read_structure(mol, input_unit, input_format, error)
   else
      call read_structure(mol, input, error, input_format)
   end if
   if (allocated(error)) then
      write(error_unit, '(a)') error%message
      error stop
   end if

   lc_method = lowercase(method)
   allocate(gradient(3, mol%nat))
   allocate(ifrez(mol%nat), source=0)

   if (echo) then
      call header(output_unit)
   endif

   if (.not.allocated(mol%lattice)) then
      lattice = 0.0_wp
   else
      lattice = transpose(mol%lattice)
   end if

   call gcp_call(mol%nat, mol%xyz, lattice, mol%num(mol%id), &
      & energy, gradient, gradlatt, dograd, dohess, any(mol%periodic), &
      & lc_method, echo, parfile)

   if (echo) then
      if (dograd) then
         call info(output_unit, energy, gradient)
      else
         call info(output_unit, energy)
      end if
   end if

    if (dograd) then
       call wregrad_tm(mol%nat, mol%nat, mol%xyz, mol%num(mol%id), ifrez, &
          & energy, gradient, echo)
    endif

contains


subroutine help(unit)
   integer, intent(in) :: unit

   write(unit, '(a, *(1x, a))') &
      "Usage: "//prog_name//" [options] <input>"

   write(unit, '(a)') &
      "", &
      ""

   write(unit, '(2x, a, t25, a)') &
      "-i, --input <format>", "Hint for the format of the input file", &
      "-l, --level <method>", "specify method to corrcet for", &
      "", "if not present parameters are read from ~/.gcppar.$HOSTNAME instead", &
      "--grad", "request gradient evaluation", &
      "--hess", "request hessian evaluation", &
      "--noprint", "Reduce printout, only print warnings", &
      "--parfile", "write gcp.param file", &
      "--version", "Print program version and exit", &
      "--help", "Show this help message"

   write(unit, '(a)')

end subroutine help


subroutine header(unit)
   integer, intent(in) :: unit

   write(unit, '(a)') &
      & " -------------------------------------------", &
      & "|              **  g C P  **                |", &
      & "|  a geometrical counterpoise correction    |", &
      & "|     H.Kruse J.G.Brandenburg S.Grimme      |", &
      & " -------------------------------------------"
   call version(unit)
   write(unit, '(a)') &
      & "", &
      & "Please cite work done with this code as:", &
      & "H. Kruse, S. Grimme J. Chem. Phys. 136, 154101 (2012)", &
      & "DOI: 10.1063/1.3700154", &
      & "For the periodic version, please also cite:", &
      & "J. G. Brandenburg, M. Alessio, B. Civalleri, M. F. Peintinger", &
      & "T. Bredow, S.Grimme J. Phys. Chem. A 117, 9282-9292 (2013).", &
      & "DOI: 10.1021/jp406658y", &
      & ""

end subroutine header


subroutine info(unit, energy, gradient)
   integer, intent(in) :: unit
   real(wp), intent(in) :: energy
   real(wp), intent(in), optional :: gradient(:, :)

   write(unit,'(a)')
   write(unit,'(a)') '** gCP correction ** '
   write(unit,'(2x,a7,F18.10,'' / (a.u.) || '',x,F11.4,'' / (kcal/mol)'')') &
      & 'Egcp:  ', energy,energy*627.5099d0
   write(unit,'(a)')
   if (present(gradient)) then
      write(unit,*)'|G|=',sum(abs(gradient))
   endif
end subroutine info


subroutine version(unit)
   integer, intent(in) :: unit
   character(len=:), allocatable :: version_string

   call get_gcp_version(string=version_string)
   write(unit, '(a, *(1x, a))') &
      & prog_name, "version", version_string

end subroutine version


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


subroutine get_arguments(input, input_format, method, dograd, dohess, echo, &
      & parfile, error)

   !> Input file name
   character(len=:), allocatable :: input

   !> Input file format
   integer, allocatable, intent(out) :: input_format

   !> Method name
   character(len=:), allocatable :: method

   !> Perform gradient calculation
   logical :: dograd

   !> Perform hessian calculation
   logical :: dohess

   !> Print information
   logical :: echo

   !> Use parameter file
   logical :: parfile

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iarg, narg
   character(len=:), allocatable :: arg

   dograd = .false.
   dohess = .false.
   echo = .true.
   parfile = .false.
   iarg = 0
   narg = command_argument_count()
   do while(iarg < narg)
      iarg = iarg + 1
      call get_argument(iarg, arg)
      select case(arg)
      case("-help", "--help")
         call help(output_unit)
         stop
      case("-version", "--version")
         call version(output_unit)
         stop
      case default
         if (.not.allocated(input)) then
            call move_alloc(arg, input)
            cycle
         end if
         call fatal_error(error, "Too many positional arguments present")
         exit
      case("-i", "-input", "--input")
         iarg = iarg + 1
         call get_argument(iarg, arg)
         if (.not.allocated(arg)) then
            call fatal_error(error, "Missing argument for input format")
            exit
         end if
         input_format = get_filetype("."//arg)
      case("-grad", "--grad")
         dograd = .true.
      case("-hess", "--hess")
         dohess = .true.
      case("-noprint", "--noprint")
         echo = .false.
      case("-parfile", "--parfile")
         parfile = .true.
      case("-l", "-level", "--level", "-func", "--func")
         iarg = iarg + 1
         call get_argument(iarg, arg)
         if (.not.allocated(arg)) then
            call fatal_error(error, "Missing argument for method")
            exit
         end if
         call move_alloc(arg, method)
      end select
   end do

   if (.not.allocated(method)) then
      method = ''
   end if

   if (.not.allocated(input)) then
      if (.not.allocated(error)) then
         call help(output_unit)
         error stop
      end if
   end if

end subroutine get_arguments


end program main
