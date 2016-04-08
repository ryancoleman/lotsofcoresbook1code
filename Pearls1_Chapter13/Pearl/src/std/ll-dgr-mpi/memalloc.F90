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
!! - \b p, \b p_f     (pressure)
!! - \b u, \b u_f     (velocity)
!! - \b rho, \b rho_f (density)
!! - \b psi_f         (chemical potential)
!! - \b ni, \b ni_f   (near neighbors)
!! - \b f, \b fbar    (order parameter distribution function)
!! - \b g, \b gbar    (momentum distribution function)
!! - \b GX, \b GY     (components of the interfacial force - potential formulation)
!! - \b gradRhoX, \b gradRhoX_f   (x component of the gradient of the density)
!! - \b gradRhoY, \b gradRhoY_f   (y component of the gradient of the density)
!! - \b gradRhoSq_f  (gradient of the density squared)
!! - \b gradRhoXX_f  (x component of the gradient of the density squared)
!! - \b gradRhoXY_f  (product of the x and y components of the density gradient)
!! - \b gradRhoYY_f  (y component of the gradient of the density squared)
!!
!! and the buffer arrays required by the parallel dual grid D2Q9 Lee-Lin
!! multiphase LBM. Arrays are deallocated when FLAG is different from one.

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
 USE FluidParams
 USE LBMParams
 USE MPIParams
 IMPLICIT NONE

! Input/Output Parameters
 INTEGER, INTENT (IN) :: FLAG

 IF(FLAG == 1) THEN
! Memory allocation for common arrays (pressure grid)
   ALLOCATE(    p(xlg:xug,ylg:yug) )
   ALLOCATE(  rho(xlg2:xug2,ylg2:yug2) )
   ALLOCATE(    u(xlg:xug,ylg:yug,1:2) )
   ALLOCATE(   ni(xlg:xug,ylg:yug,1:4) )
   ALLOCATE(    g(xlg:xug,ylg:yug,0:gdim) )
   ALLOCATE( gbar(xlg:xug,ylg:yug,0:gdim) )

! Memory allocation for common arrays (order parameter grid)
   ALLOCATE(   p_f(xlg_f:xug_f,ylg_f:yug_f) )
   ALLOCATE( rho_f(xlg3_f:xug3_f,ylg3_f:yug3_f) )
   ALLOCATE( psi_f(xlg2_f:xug2_f,ylg2_f:yug2_f) )
   ALLOCATE(   u_f(xlg_f:xug_f,ylg_f:yug_f,1:2) )
   ALLOCATE(  ni_f(xlg2_f:xug2_f,ylg2_f:yug2_f,1:4) )
   ALLOCATE(     f(xlg_f:xug_f,ylg_f:yug_f,0:fdim) )
   ALLOCATE(  fbar(xlg_f:xug_f,ylg_f:yug_f,0:fdim) )

! Memory allocation for differential terms (pressure grid)
   ALLOCATE(       GX(xlg:xug,ylg:yug) )
   ALLOCATE(       GY(xlg:xug,ylg:yug) )
   ALLOCATE( gradRhoX(xlg:xug,ylg:yug) )
   ALLOCATE( gradRhoY(xlg:xug,ylg:yug) )

! Memory allocation for differential terms (momentum calculation)
   ALLOCATE(  gradRhoX_f(xlg_f:xug_f,ylg_f:yug_f) )
   ALLOCATE(  gradRhoY_f(xlg_f:xug_f,ylg_f:yug_f) )
   ALLOCATE( gradRhoSq_f(xlg_f:xug_f,ylg_f:yug_f) )
   ALLOCATE( gradRhoXX_f(xlg_f:xug_f,ylg_f:yug_f) )
   ALLOCATE( gradRhoXY_f(xlg_f:xug_f,ylg_f:yug_f) )
   ALLOCATE( gradRhoYY_f(xlg_f:xug_f,ylg_f:yug_f) )


! Allocation for the exchanged data in the x direction
   xsize    = yu   - yl   + 1
   xsize_f  = yu_f - yl_f + 1
   xsize2   = xsize   + xsize
   xsize3   = xsize   + xsize   + xsize
   xsize3_f = xsize_f + xsize_f + xsize_f

   ALLOCATE( rho_west_snd(1:xsize3_f) )
   ALLOCATE( rho_east_snd(1:xsize3_f) )
   ALLOCATE( rho_west_rcv(1:xsize3_f) )
   ALLOCATE( rho_east_rcv(1:xsize3_f) )

   ALLOCATE( rhoG_west_snd(1:xsize2) )
   ALLOCATE( rhoG_east_snd(1:xsize2) )
   ALLOCATE( rhoG_west_rcv(1:xsize2) )
   ALLOCATE( rhoG_east_rcv(1:xsize2) )

   ALLOCATE( f_west_snd(1:xsize3_f) )
   ALLOCATE( f_east_snd(1:xsize3_f) )
   ALLOCATE( f_west_rcv(1:xsize3_f) )
   ALLOCATE( f_east_rcv(1:xsize3_f) )

   ALLOCATE( hydro_west_snd(1:xsize3) )
   ALLOCATE( hydro_east_snd(1:xsize3) )
   ALLOCATE( hydro_west_rcv(1:xsize3) )
   ALLOCATE( hydro_east_rcv(1:xsize3) )

   ALLOCATE( g_west_snd(1:xsize3) )
   ALLOCATE( g_east_snd(1:xsize3) )
   ALLOCATE( g_west_rcv(1:xsize3) )
   ALLOCATE( g_east_rcv(1:xsize3) )

! Allocation for the exchanged data in the y direction
   ysize    = xu   - xl   + 1
   ysize_f  = xu_f - xl_f + 1
   ysize2   = ysize   + ysize
   ysize3   = ysize   + ysize   + ysize
   ysize3_f = ysize_f + ysize_f + ysize_f
   
   ALLOCATE( rho_south_snd(1:ysize3_f) )
   ALLOCATE( rho_north_snd(1:ysize3_f) )
   ALLOCATE( rho_south_rcv(1:ysize3_f) )
   ALLOCATE( rho_north_rcv(1:ysize3_f) )

   ALLOCATE( rhoG_south_snd(1:ysize2) )
   ALLOCATE( rhoG_north_snd(1:ysize2) )
   ALLOCATE( rhoG_south_rcv(1:ysize2) )
   ALLOCATE( rhoG_north_rcv(1:ysize2) )

   ALLOCATE( f_south_snd(1:ysize3_f) )
   ALLOCATE( f_north_snd(1:ysize3_f) )
   ALLOCATE( f_south_rcv(1:ysize3_f) )
   ALLOCATE( f_north_rcv(1:ysize3_f) )

   ALLOCATE( hydro_south_snd(1:ysize3) )
   ALLOCATE( hydro_north_snd(1:ysize3) )
   ALLOCATE( hydro_south_rcv(1:ysize3) )
   ALLOCATE( hydro_north_rcv(1:ysize3) )

   ALLOCATE( g_south_snd(1:ysize3) )
   ALLOCATE( g_north_snd(1:ysize3) )
   ALLOCATE( g_south_rcv(1:ysize3) )
   ALLOCATE( g_north_rcv(1:ysize3) )

 ELSE
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
   DEALLOCATE( hydro_west_snd )
   DEALLOCATE( hydro_west_rcv )
   DEALLOCATE( hydro_east_snd )
   DEALLOCATE( hydro_east_rcv )
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
   DEALLOCATE( hydro_south_snd )
   DEALLOCATE( hydro_south_rcv )
   DEALLOCATE( hydro_north_snd )
   DEALLOCATE( hydro_north_rcv )

   DEALLOCATE( p )
   DEALLOCATE( u )
   DEALLOCATE( g )
   DEALLOCATE( ni )
   DEALLOCATE( rho )
   DEALLOCATE( gbar )
   DEALLOCATE( gradRhoX )
   DEALLOCATE( gradRhoY )
   DEALLOCATE( GX )
   DEALLOCATE( GY )

   DEALLOCATE( p_f )
   DEALLOCATE( u_f )
   DEALLOCATE( f )
   DEALLOCATE( ni_f )
   DEALLOCATE( rho_f )
   DEALLOCATE( psi_f )
   DEALLOCATE( fbar )
   DEALLOCATE( gradRhoX_f )
   DEALLOCATE( gradRhoY_f )
   DEALLOCATE( gradRhoSq_f )
   DEALLOCATE( gradRhoXX_f )
   DEALLOCATE( gradRhoXY_f )
   DEALLOCATE( gradRhoYY_f )
 END IF

 RETURN
 END SUBROUTINE MemAlloc
