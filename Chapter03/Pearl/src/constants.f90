!Copyright (c) 2014, Per Berg and Jacob Weismann Poulsen, DMI
!All rights reserved.
!
!Redistribution and use in source and binary forms, with or without
!modification, are permitted provided that the following conditions are met: 
! 
!1. Redistributions of source code must retain the above copyright notice, this
!   list of conditions and the following disclaimer. 
!2. Redistributions in binary form must reproduce the above copyright notice,
!   this list of conditions and the following disclaimer in the documentation
!   and/or other materials provided with the distribution. 
!
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
!ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
!WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
!ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
!(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
!ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
!(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
! 
!The views and conclusions contained in the software and documentation are those
!of the authors and should not be interpreted as representing official policies,
!either expressed or implied, of the FreeBSD Project.

module constants

!-------------------------------------------------------------------------------
 
  !- directives ----------------------------------------------------------------
  implicit none
  private

  !- constants -----------------------------------------------------------------
  real(8), parameter, public  :: g      =    9.81_8
  real(8), parameter, public  :: qrt = 0.25_8, tenth = 0.1_8
  real(8), parameter, public  :: hgr    =    0.05_8
  real(8), parameter, public  :: hnull  =    0.1_8
  real(8), parameter, public  :: heps   = 0.00001_8
  !-----------------------------------------------------------------------------


  !- time step information -----------------------------------------------------
  integer(4),    save, public :: izeits, jahra
  real(8),       save, public :: rzeite, zzeits, tim, zeitstu
  character(11), save, public :: aufdat, starta
  !-----------------------------------------------------------------------------

  !- optional parameter list ---------------------------------------------------
  real(8),    save, public :: decomp_coeff(3)
  integer(4), save, public :: nclss
  integer(4), save, public :: itmax, allowsubdiv
  logical,    save, public :: ionml
  integer(4), save, public :: iot=10000 
  logical,    save, public :: ompsplit_exact=.false., ompsplit_op=.false.
  logical,    save, public :: ompsplit_load=.false., omp_debug=.true.
  logical,    save, public :: only_islices, omp_loadbalance_experiment=.false.
  logical,    save, public :: coldst, spinup, archive, ldebug
  real(8), parameter, public  :: zero = 0.0_8, one = 1.0_8, half = 0.5_8
  real(8), parameter, public  :: two = 2.0_8, three = 3.0_8, four = 4.0_8
  real(8), parameter, public  :: onethird=1.0_8/3.0_8, twothird=2.0_8/3.0_8
  real(8), parameter, public  :: fourthird=4.0_8/3.0_8
  real(8), parameter, public  :: cpw    = 4180.00_8
  real(8), parameter, public  :: rhow   = 1027.00_8
  logical,    save, public :: autodecomp_op=.false.
  integer(4), save, public :: nproci=0, nprocj=0, decomp_version=1
  logical,    save, public :: msscrr=.false.
end module constants
