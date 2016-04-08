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
!! - \b p          (pressure)
!! - \b u          (velocity)
!! - \b rho        (density)
!! - \b psi        (chemical potential)
!! - \b ni         (near neighbors)
!! - \b f, \b fbar (order parameter distribution function)
!! - \b g, \b gbar (momentum distribution function)
!! - \b gradRhoX   (x component of the gradient of the density)
!! - \b gradRhoY   (y component of the gradient of the density)
!! - \b gradRhoSq  (gradient of the density squared)
!! - \b gradRhoXX  (x component of the gradient of the density squared)
!! - \b gradRhoXY  (product of the x and y components of the density gradient)
!! - \b gradRhoYY  (y component of the gradient of the density squared)
!!
!! and the buffer arrays required by the parallel D2Q9 Lee-Lin multiphase LBM.
!! Arrays are deallocated when FLAG is different from one.

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
 USE FluidParams
 USE LBMParams,  ONLY : f, fbar, fdim, g, gbar, gdim
 USE MPIParams
 IMPLICIT NONE

 !   Input/Output Parameters
 INTEGER, INTENT (IN) :: FLAG

 IF(FLAG == 1) THEN

! Allocation for hydrodynamics arrays
   ALLOCATE(    p(xlg:xug,ylg:yug) )
   ALLOCATE(  rho(xlg3:xug3,ylg3:yug3) )
   ALLOCATE(  psi(xlg2:xug2,ylg2:yug2) )
   ALLOCATE(    u(xlg:xug,ylg:yug,1:2) )
   ALLOCATE(   ni(xlg2:xug2,ylg2:yug2,1:4) )

! Allocation for distribution functions
   ALLOCATE(    f(xlg:xug,ylg:yug,0:fdim) )
   ALLOCATE(    g(xlg:xug,ylg:yug,0:gdim) )
   ALLOCATE( fbar(xlg:xug,ylg:yug,0:fdim) )
   ALLOCATE( gbar(xlg:xug,ylg:yug,0:gdim) )

! Allocation for differential terms
   ALLOCATE(  gradRhoX(xlg:xug,ylg:yug) )
   ALLOCATE(  gradRhoY(xlg:xug,ylg:yug) )
   ALLOCATE( gradRhoSq(xlg:xug,ylg:yug) )
   ALLOCATE( gradRhoXX(xlg:xug,ylg:yug) )
   ALLOCATE( gradRhoXY(xlg:xug,ylg:yug) )
   ALLOCATE( gradRhoYY(xlg:xug,ylg:yug) )

! Allocation for the exchanged data in the x direction
   xsize  = yu - yl + 1
   xsize2 = xsize + xsize
   xsize3 = xsize2 + xsize

   ALLOCATE( rho_west_snd(1:xsize3) )
   ALLOCATE( rho_east_snd(1:xsize3) )
   ALLOCATE( rho_west_rcv(1:xsize3) )
   ALLOCATE( rho_east_rcv(1:xsize3) )

   ALLOCATE( f_west_snd(1:xsize3) )
   ALLOCATE( f_east_snd(1:xsize3) )
   ALLOCATE( f_west_rcv(1:xsize3) )
   ALLOCATE( f_east_rcv(1:xsize3) )

   ALLOCATE( g_west_snd(1:xsize3) )
   ALLOCATE( g_east_snd(1:xsize3) )
   ALLOCATE( g_west_rcv(1:xsize3) )
   ALLOCATE( g_east_rcv(1:xsize3) )

! Allocation for the exchanged data in the y direction
   ysize  = xu - xl + 1
   ysize2 = ysize + ysize
   ysize3 = ysize2 + ysize

   ALLOCATE( rho_south_snd(1:ysize3) )
   ALLOCATE( rho_north_snd(1:ysize3) )
   ALLOCATE( rho_south_rcv(1:ysize3) )
   ALLOCATE( rho_north_rcv(1:ysize3) )

   ALLOCATE( f_south_snd(1:ysize3) )
   ALLOCATE( f_north_snd(1:ysize3) )
   ALLOCATE( f_south_rcv(1:ysize3) )
   ALLOCATE( f_north_rcv(1:ysize3) )

   ALLOCATE( g_south_snd(1:ysize3) )
   ALLOCATE( g_north_snd(1:ysize3) )
   ALLOCATE( g_south_rcv(1:ysize3) )
   ALLOCATE( g_north_rcv(1:ysize3) )

 ELSE

   DEALLOCATE( p )
   DEALLOCATE( u )
   DEALLOCATE( ni )
   DEALLOCATE( rho )
   DEALLOCATE( psi )

   DEALLOCATE( f )
   DEALLOCATE( g )
   DEALLOCATE( fbar )
   DEALLOCATE( gbar )

   DEALLOCATE( gradRhoX )
   DEALLOCATE( gradRhoY )
   DEALLOCATE( gradRhoSq )
   DEALLOCATE( gradRhoXX )
   DEALLOCATE( gradRhoXY )
   DEALLOCATE( gradRhoYY )

   DEALLOCATE( f_west_snd )
   DEALLOCATE( f_west_rcv )
   DEALLOCATE( f_east_snd )
   DEALLOCATE( f_east_rcv )
   DEALLOCATE( g_west_snd )
   DEALLOCATE( g_east_snd )
   DEALLOCATE( g_west_rcv )
   DEALLOCATE( g_east_rcv )
   DEALLOCATE( rho_west_snd )
   DEALLOCATE( rho_west_rcv )
   DEALLOCATE( rho_east_snd )
   DEALLOCATE( rho_east_rcv )
   DEALLOCATE( f_south_snd )
   DEALLOCATE( f_south_rcv )
   DEALLOCATE( f_north_snd )
   DEALLOCATE( f_north_rcv )
   DEALLOCATE( g_south_snd )
   DEALLOCATE( g_south_rcv )
   DEALLOCATE( g_north_snd )
   DEALLOCATE( g_north_rcv )
   DEALLOCATE( rho_south_snd )
   DEALLOCATE( rho_south_rcv )
   DEALLOCATE( rho_north_snd )
   DEALLOCATE( rho_north_rcv )
 END IF

 RETURN
 END SUBROUTINE MemAlloc

