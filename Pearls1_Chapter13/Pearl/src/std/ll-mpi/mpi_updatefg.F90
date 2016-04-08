!-------------------------------------------------------------------------------
! Subroutine : MPI_UpdateFG
! Revision   : 1.0 (2008-06-15)
! Author     : Carlos Rosales Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
!> @file
!! Update f and g in nodes at processor boundaries using MPI and update rho.
!> @details
!! Copy outward values of the order parameter (f) and the pressure (g)
!! distribution functions from the ghosts into inward values in the real
!! nodes. The density (rho) is calculated using the updated values of the order
!! parameter distribution function f.
!!
!! MPI exchange is done using non-blocking receive and blocking send pairs
!! (completed by a wait statement) in alternate directions to avoid deadlock.
!! The components of g and f that need to be exchanged are packed into buffer
!! arrays before the exchange.

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

 SUBROUTINE MPI_UpdateFG

! Common Variables
 USE Domain,      ONLY : xl, xlg, xsize3, xu, xug, yl, ylg, ysize3, yu, yug
 USE LBMParams,   ONLY : fbar, gbar
 USE FluidParams, ONLY : rho
 USE MPIParams
 USE MPI
 IMPLICIT NONE

! Local Variables
 INTEGER :: count, i, j
 INTEGER :: MPI_ERR, MPI_REQ
 INTEGER :: status(MPI_STATUS_SIZE)

!--------- Save data for inward g components -----------------------------------
! X DIRECTION
 count = 1
 DO j = yl, yu
   g_west_snd(count)     = gbar(xlg,j,3)
   g_west_snd(count + 1) = gbar(xlg,j,6)
   g_west_snd(count + 2) = gbar(xlg,j,7)

   f_west_snd(count)     = fbar(xlg,j,3)
   f_west_snd(count + 1) = fbar(xlg,j,6)
   f_west_snd(count + 2) = fbar(xlg,j,7)

   g_east_snd(count)     = gbar(xug,j,1)
   g_east_snd(count + 1) = gbar(xug,j,5)
   g_east_snd(count + 2) = gbar(xug,j,8)

   f_east_snd(count)     = fbar(xug,j,1)
   f_east_snd(count + 1) = fbar(xug,j,5)
   f_east_snd(count + 2) = fbar(xug,j,8)
   count = count + 3
 END DO

! Y DIRECTION
 count = 1
 DO i = xl, xu
   g_south_snd(count)     = gbar(i,ylg,4)
   g_south_snd(count + 1) = gbar(i,ylg,7)
   g_south_snd(count + 2) = gbar(i,ylg,8)

   f_south_snd(count)     = fbar(i,ylg,4)
   f_south_snd(count + 1) = fbar(i,ylg,7)
   f_south_snd(count + 2) = fbar(i,ylg,8)

   g_north_snd(count)     = gbar(i,yug,2)
   g_north_snd(count + 1) = gbar(i,yug,5)
   g_north_snd(count + 2) = gbar(i,yug,6)

   f_north_snd(count)     = fbar(i,yug,2)
   f_north_snd(count + 1) = fbar(i,yug,5)
   f_north_snd(count + 2) = fbar(i,yug,6)
   count = count + 3
 END DO
 
! DIAGONALS
 g_ne_snd = gbar(xug,yug,5)
 g_nw_snd = gbar(xlg,yug,6)
 g_sw_snd = gbar(xlg,ylg,7)
 g_se_snd = gbar(xug,ylg,8)
 f_ne_snd = fbar(xug,yug,5)
 f_nw_snd = fbar(xlg,yug,6)
 f_sw_snd = fbar(xlg,ylg,7)
 f_se_snd = fbar(xug,ylg,8)

!---------Exchange information using MPI calls----------------------------------
! X DIRECTION
 CALL MPI_IRECV(f_east_rcv, xsize3, MPI_DOUBLE_PRECISION, east, TAG1, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( f_west_snd, xsize3, MPI_DOUBLE_PRECISION, west, TAG1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(g_east_rcv, xsize3, MPI_DOUBLE_PRECISION, east, TAG2, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( g_west_snd, xsize3, MPI_DOUBLE_PRECISION, west, TAG2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(f_west_rcv, xsize3, MPI_DOUBLE_PRECISION, west, TAG3, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( f_east_snd, xsize3, MPI_DOUBLE_PRECISION, east, TAG3, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(g_west_rcv, xsize3, MPI_DOUBLE_PRECISION, west, TAG4, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( g_east_snd, xsize3, MPI_DOUBLE_PRECISION, east, TAG4, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

! Y DIRECTION
 CALL MPI_IRECV(f_north_rcv, ysize3, MPI_DOUBLE_PRECISION, north, TAG1, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( f_south_snd, ysize3, MPI_DOUBLE_PRECISION, south, TAG1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(g_north_rcv, ysize3, MPI_DOUBLE_PRECISION, north, TAG2, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( g_south_snd, ysize3, MPI_DOUBLE_PRECISION, south, TAG2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(f_south_rcv, ysize3, MPI_DOUBLE_PRECISION, south, TAG3, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( f_north_snd, ysize3, MPI_DOUBLE_PRECISION, north, TAG3, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(g_south_rcv, ysize3, MPI_DOUBLE_PRECISION, south, TAG4, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( g_north_snd, ysize3, MPI_DOUBLE_PRECISION, north, TAG4, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

! DIAGONALS
 CALL MPI_IRECV(g_ne_rcv, 1, MPI_DOUBLE_PRECISION, ne, TAG1, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( g_sw_snd, 1, MPI_DOUBLE_PRECISION, sw, TAG1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(g_nw_rcv, 1, MPI_DOUBLE_PRECISION, nw, TAG2, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( g_se_snd, 1, MPI_DOUBLE_PRECISION, se, TAG2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(g_se_rcv, 1, MPI_DOUBLE_PRECISION, se, TAG3, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( g_nw_snd, 1, MPI_DOUBLE_PRECISION, nw, TAG3, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(g_sw_rcv, 1, MPI_DOUBLE_PRECISION, sw, TAG4, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( g_ne_snd, 1, MPI_DOUBLE_PRECISION, ne, TAG4, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)
 
 CALL MPI_IRECV(f_ne_rcv, 1, MPI_DOUBLE_PRECISION, ne, TAG1, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( f_sw_snd, 1, MPI_DOUBLE_PRECISION, sw, TAG1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(f_nw_rcv, 1, MPI_DOUBLE_PRECISION, nw, TAG2, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( f_se_snd, 1, MPI_DOUBLE_PRECISION, se, TAG2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(f_se_rcv, 1, MPI_DOUBLE_PRECISION, se, TAG3, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( f_nw_snd, 1, MPI_DOUBLE_PRECISION, nw, TAG3, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(f_sw_rcv, 1, MPI_DOUBLE_PRECISION, sw, TAG4, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( f_ne_snd, 1, MPI_DOUBLE_PRECISION, ne, TAG4, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)


!---------Update inward components of f and g in the real nodes-----------------
! X DIRECTION
 count = 1
 DO j = yl, yu
   gbar(xu,j,3) = g_east_rcv(count)
   gbar(xu,j,6) = g_east_rcv(count + 1)
   gbar(xu,j,7) = g_east_rcv(count + 2)

   fbar(xu,j,3) = f_east_rcv(count)
   fbar(xu,j,6) = f_east_rcv(count + 1)
   fbar(xu,j,7) = f_east_rcv(count + 2)

   gbar(xl,j,1) = g_west_rcv(count)
   gbar(xl,j,5) = g_west_rcv(count + 1)
   gbar(xl,j,8) = g_west_rcv(count + 2)

   fbar(xl,j,1) = f_west_rcv(count)
   fbar(xl,j,5) = f_west_rcv(count + 1)
   fbar(xl,j,8) = f_west_rcv(count + 2)
   count = count + 3
 END DO

! Y DIRECTION
 count = 1
 DO i = xl, xu
   gbar(i,yu,4) = g_north_rcv(count)
   gbar(i,yu,7) = g_north_rcv(count + 1)
   gbar(i,yu,8) = g_north_rcv(count + 2)

   fbar(i,yu,4) = f_north_rcv(count)
   fbar(i,yu,7) = f_north_rcv(count + 1)
   fbar(i,yu,8) = f_north_rcv(count + 2)

   gbar(i,yl,2) = g_south_rcv(count)
   gbar(i,yl,5) = g_south_rcv(count + 1)
   gbar(i,yl,6) = g_south_rcv(count + 2)

   fbar(i,yl,2) = f_south_rcv(count)
   fbar(i,yl,5) = f_south_rcv(count + 1)
   fbar(i,yl,6) = f_south_rcv(count + 2)
   count = count + 3
 END DO

! DIAGONALS
 gbar(xu,yu,7) = g_ne_rcv
 gbar(xl,yu,8) = g_nw_rcv
 gbar(xl,yl,5) = g_sw_rcv
 gbar(xu,yl,6) = g_se_rcv
 fbar(xu,yu,7) = f_ne_rcv
 fbar(xl,yu,8) = f_nw_rcv
 fbar(xl,yl,5) = f_sw_rcv
 fbar(xu,yl,6) = f_se_rcv


!--------- Calculate the updated density ---------------------------------------
 DO j = yl, yu
   DO i = xl, xu

     rho(i,j) = fbar(i,j,0) + fbar(i,j,1) + fbar(i,j,2) + fbar(i,j,3) &
              + fbar(i,j,4) + fbar(i,j,5) + fbar(i,j,6) + fbar(i,j,7) &
              + fbar(i,j,8)

   END DO
 END DO


 RETURN
 END SUBROUTINE MPI_UpdateFG

