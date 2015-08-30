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
!! - \b p, \b p_f   (pressure)
!! - \b u, \b u_f   (velocity)
!! - \b phi_f       (order parameter)
!! - \b ni, \b ni_f (near neighbors)
!! - \b f, \b fcol  (order parameter distribution function)
!! - \b g           (momentum distribution function)
!! - \b gradPhiX_f  (x component of order parameter gradient)
!! - \b gradPhiY_f  (y component of order parameter gradient)
!! - \b lapPhi_f    (Laplacian of the order parameter)
!!
!! and the MPI buffer arrays required by the dual grid parallel D2Q5/D2Q9
!! Zheng-Shu-Chew multiphase LBM. Memory is deallocated when FLAG is different
!! from one.

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

! Common Variables
 USE Domain
 USE FluidParams, ONLY : gradPhiX_f, gradPhiY_f, lapPhi_f, p, p_f, phi_f, u, u_f
 USE LBMParams,   ONLY : f, fcol, fdim, g, gdim
 USE MPIParams
 IMPLICIT NONE

! Input/Output Parameters
 INTEGER, INTENT (IN) :: FLAG

!-------- Allocate memory for all the MPI data exchanges -----------------------
 IF(FLAG == 1) THEN

! Memory allocation - Momentum/Pressure Mesh
   ALLOCATE(  p(xlg:xug,ylg:yug) )
   ALLOCATE(  u(xlg:xug,ylg:yug,2) )
   ALLOCATE(  g(xlg:xug,ylg:yug,0:gdim,0:1) )
   ALLOCATE( ni(xlg:xug,ylg:yug,1:4) )

! Memory Allocation - Order Parameter Mesh
   ALLOCATE(   p_f(xlg_f:xug_f,ylg_f:yug_f) )
   ALLOCATE(   u_f(xlg_f:xug_f,ylg_f:yug_f,2) )
   ALLOCATE( phi_f(xlg_f:xug_f,ylg_f:yug_f) )
   ALLOCATE(     f(xlg_f:xug_f,ylg_f:yug_f,0:fdim) )
   ALLOCATE(  fcol(xlg_f:xug_f,ylg_f:yug_f,0:fdim) )
   ALLOCATE(  ni_f(xlg_f:xug_f,ylg_f:yug_f,1:4) )
   ALLOCATE( gradPhiX_f(xlg_f:xug_f,ylg_f:yug_f) )
   ALLOCATE( gradPhiY_f(xlg_f:xug_f,ylg_f:yug_f) )
   ALLOCATE(   lapPhi_f(xlg_f:xug_f,ylg_f:yug_f) )

! Allocation for the exchanged data in the x direction
   xsize_f = yu_f - yl_f + 1
   xsize   = yu - yl + 1
   xsize6  = xsize*6

   ALLOCATE( phi_west_snd(1:xsize_f) )
   ALLOCATE( phi_east_snd(1:xsize_f) )
   ALLOCATE( phi_west_rcv(1:xsize_f) )
   ALLOCATE( phi_east_rcv(1:xsize_f) )

   ALLOCATE( f_west_snd(1:xsize_f) )
   ALLOCATE( f_east_snd(1:xsize_f) )
   ALLOCATE( f_west_rcv(1:xsize_f) )
   ALLOCATE( f_east_rcv(1:xsize_f) )

   ALLOCATE( hydro_west_snd(1:xsize6) )
   ALLOCATE( hydro_east_snd(1:xsize6) )
   ALLOCATE( hydro_west_rcv(1:xsize6) )
   ALLOCATE( hydro_east_rcv(1:xsize6) )

! Allocation for the exchanged data in the y direction
   ysize_f = xu_f - xl_f + 1
   ysize   = xu - xl + 1
   ysize6  = ysize*6

   ALLOCATE( phi_south_snd(1:ysize_f) )
   ALLOCATE( phi_north_snd(1:ysize_f) )
   ALLOCATE( phi_south_rcv(1:ysize_f) )
   ALLOCATE( phi_north_rcv(1:ysize_f) )

   ALLOCATE( f_south_snd(1:ysize_f) )
   ALLOCATE( f_north_snd(1:ysize_f) )
   ALLOCATE( f_south_rcv(1:ysize_f) )
   ALLOCATE( f_north_rcv(1:ysize_f) )

   ALLOCATE( hydro_south_snd(1:ysize6) )
   ALLOCATE( hydro_north_snd(1:ysize6) )
   ALLOCATE( hydro_south_rcv(1:ysize6) )
   ALLOCATE( hydro_north_rcv(1:ysize6) )

!-------- Deallocate memory corresponding to MPI data exchanges ----------------
 ELSE
   DEALLOCATE( f )
   DEALLOCATE( fcol )
   DEALLOCATE( g )
   DEALLOCATE( ni )
   DEALLOCATE( p )
   DEALLOCATE( u )
   DEALLOCATE( ni_f )
   DEALLOCATE( phi_f )
   DEALLOCATE( p_f )
   DEALLOCATE( u_f )
   DEALLOCATE( gradPhiX_f )
   DEALLOCATE( gradPhiY_f )
   DEALLOCATE( lapPhi_f )
   DEALLOCATE( f_west_snd )
   DEALLOCATE( f_west_rcv )
   DEALLOCATE( f_east_snd )
   DEALLOCATE( f_east_rcv )
   DEALLOCATE( phi_west_snd )
   DEALLOCATE( phi_west_rcv )
   DEALLOCATE( phi_east_snd )
   DEALLOCATE( phi_east_rcv )
   DEALLOCATE( hydro_west_snd )
   DEALLOCATE( hydro_east_snd )
   DEALLOCATE( hydro_west_rcv )
   DEALLOCATE( hydro_east_rcv )
   DEALLOCATE( f_south_snd )
   DEALLOCATE( f_south_rcv )
   DEALLOCATE( f_north_snd )
   DEALLOCATE( f_north_rcv )
   DEALLOCATE( phi_south_snd )
   DEALLOCATE( phi_south_rcv )
   DEALLOCATE( phi_north_snd )
   DEALLOCATE( phi_north_rcv )
   DEALLOCATE( hydro_south_snd )
   DEALLOCATE( hydro_south_rcv )
   DEALLOCATE( hydro_north_snd )
   DEALLOCATE( hydro_north_rcv )
 END IF

 RETURN
 END SUBROUTINE MemAlloc
