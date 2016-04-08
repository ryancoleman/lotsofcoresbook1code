!-------------------------------------------------------------------------------
! Subroutine : MemAlloc
! Revision   : 1.0 (2008-06-15)
! Author     : Carlos Rosales Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
!> @file
!! Allocate (FLAG == 1) or deallocate (FLAG != 1) memory for common arrays.
!> @details
!! If FLAG is 1 then memory is allocated for the following common arrays:
!!
!! - \b p   (pressure)
!! - \b u   (velocity)
!! - \b phi (order parameter)
!! - \b ni  (near neighbors)
!! - \b f   (order parameter distribution function)
!! - \b g   (momentum distribution function)
!!
!! and the MPI buffer arrays required by the parallel D2Q5/D2Q9 Zheng-Shu-Chew
!! multiphase LBM. Memory is deallocated when FLAG is different from one.

!-------------------------------------------------------------------------------
! Copyright 2008 Carlos Rosales Fernandez, David S. Whyte and IHPC (A*STAR).
!
! This file is part of MP-LABS.
!
! MP-LABS is free software: you can redistribute it and/or modify it under the
! terms of the GNU GPL version 3 or (at your option) any later version.
!
! MP-LABS is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
! A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with
! MP-LABS, in the file COPYING.txt. If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------

 SUBROUTINE MemAlloc(FLAG)

!  Common Variables
 USE Domain
 USE FluidParams, ONLY : p, phi, u
 USE LBMParams,   ONLY : f, fdim, g, gdim
 USE MPIParams
 IMPLICIT NONE

 !   Input/Output Parameters
 INTEGER, INTENT (IN) :: FLAG


 IF(FLAG == 1) THEN

! Hydrodynamics arrays
   ALLOCATE( p(xl:xu,yl:yu) )
   ALLOCATE( u(xl:xu,yl:yu,2) )
   ALLOCATE( phi(xlg:xug,ylg:yug) )

! Distribution function arrays
   ALLOCATE( f(xlg:xug,ylg:yug,0:fdim,0:1) )
   ALLOCATE( g(xlg:xug,ylg:yug,0:gdim,0:1) )

! Neighbor list array
   ALLOCATE( ni(xlg:xug,ylg:yug,1:4) )

! Allocation for exchange of data in the x direction
   xsize  = yu - yl + 1
   xsize2 = xsize + xsize
   xsize3 = xsize2 + xsize

   ALLOCATE( phi_west_snd(1:xsize) )
   ALLOCATE( phi_east_snd(1:xsize) )
   ALLOCATE( phi_west_rcv(1:xsize) )
   ALLOCATE( phi_east_rcv(1:xsize) )

   ALLOCATE( f_west_snd(1:xsize) )
   ALLOCATE( f_east_snd(1:xsize) )
   ALLOCATE( f_west_rcv(1:xsize) )
   ALLOCATE( f_east_rcv(1:xsize) )

   ALLOCATE( g_west_snd(1:xsize3) )
   ALLOCATE( g_east_snd(1:xsize3) )
   ALLOCATE( g_west_rcv(1:xsize3) )
   ALLOCATE( g_east_rcv(1:xsize3) )

! Allocation for exchange of data in the y direction
   ysize  = xu - xl + 1
   ysize2 = ysize + ysize
   ysize3 = ysize2 + ysize

   ALLOCATE( phi_south_snd(1:ysize) )
   ALLOCATE( phi_north_snd(1:ysize) )
   ALLOCATE( phi_south_rcv(1:ysize) )
   ALLOCATE( phi_north_rcv(1:ysize) )

   ALLOCATE( f_south_snd(1:ysize) )
   ALLOCATE( f_north_snd(1:ysize) )
   ALLOCATE( f_south_rcv(1:ysize) )
   ALLOCATE( f_north_rcv(1:ysize) )

   ALLOCATE( g_south_snd(1:ysize3) )
   ALLOCATE( g_north_snd(1:ysize3) )
   ALLOCATE( g_south_rcv(1:ysize3) )
   ALLOCATE( g_north_rcv(1:ysize3) )

 ELSE

   DEALLOCATE( f )
   DEALLOCATE( g )
   DEALLOCATE( ni )
   DEALLOCATE( phi )
   DEALLOCATE( p )
   DEALLOCATE( u )
   DEALLOCATE( f_west_snd )
   DEALLOCATE( f_west_rcv )
   DEALLOCATE( f_east_snd )
   DEALLOCATE( f_east_rcv )
   DEALLOCATE( g_west_snd )
   DEALLOCATE( g_east_snd )
   DEALLOCATE( g_west_rcv )
   DEALLOCATE( g_east_rcv )
   DEALLOCATE( phi_west_snd )
   DEALLOCATE( phi_west_rcv )
   DEALLOCATE( phi_east_snd )
   DEALLOCATE( phi_east_rcv )
   DEALLOCATE( f_south_snd )
   DEALLOCATE( f_south_rcv )
   DEALLOCATE( f_north_snd )
   DEALLOCATE( f_north_rcv )
   DEALLOCATE( g_south_snd )
   DEALLOCATE( g_south_rcv )
   DEALLOCATE( g_north_snd )
   DEALLOCATE( g_north_rcv )
   DEALLOCATE( phi_south_snd )
   DEALLOCATE( phi_south_rcv )
   DEALLOCATE( phi_north_snd )
   DEALLOCATE( phi_north_rcv )

 END IF

 RETURN
 END SUBROUTINE MemAlloc
