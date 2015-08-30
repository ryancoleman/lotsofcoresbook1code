!-------------------------------------------------------------------------------
! Subroutine : MPI_UpdateG
! Revision   : 1.0 (2008-06-15)
! Author     : Carlos Rosales Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
!> @file
!! Update gbar in nodes at processor boundaries using MPI.
!> @details
!! Copy outward values of the pressure distribution function gbar from the
!! ghosts into inward values in the real nodes after the streaming step.
!!
!! MPI exchange is done using non-blocking receive and blocking send pairs
!! (completed by a wait statement) in alternate directions to avoid deadlock.
!! The components of gbar that need to be exchanged are packed into buffer
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

 SUBROUTINE MPI_UpdateG

! Common Variables
 USE Domain
 USE FluidParams
 USE LBMParams
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

   g_east_snd(count)     = gbar(xug,j,1)
   g_east_snd(count + 1) = gbar(xug,j,5)
   g_east_snd(count + 2) = gbar(xug,j,8)
   count = count + 3
 END DO

! Y DIRECTION
 count = 1
 DO i = xl, xu
   g_south_snd(count)     = gbar(i,ylg,4)
   g_south_snd(count + 1) = gbar(i,ylg,7)
   g_south_snd(count + 2) = gbar(i,ylg,8)

   g_north_snd(count)     = gbar(i,yug,2)
   g_north_snd(count + 1) = gbar(i,yug,5)
   g_north_snd(count + 2) = gbar(i,yug,6)
   count = count + 3
 END DO
 
! DIAGONALS
 g_ne_snd = gbar(xug,yug,5)
 g_nw_snd = gbar(xlg,yug,6)
 g_sw_snd = gbar(xlg,ylg,7)
 g_se_snd = gbar(xug,ylg,8)


!--------- Exchange information using MPI calls --------------------------------
! X DIRECTION
 CALL MPI_IRECV(g_east_rcv, xsize3, MPI_DOUBLE_PRECISION, east, tag2, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( g_west_snd, xsize3, MPI_DOUBLE_PRECISION, west, tag2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(g_west_rcv, xsize3, MPI_DOUBLE_PRECISION, west, tag4, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( g_east_snd, xsize3, MPI_DOUBLE_PRECISION, east, tag4, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

! Y DIRECTION
 CALL MPI_IRECV(g_north_rcv, ysize3, MPI_DOUBLE_PRECISION, north, tag2, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( g_south_snd, ysize3, MPI_DOUBLE_PRECISION, south, tag2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(g_south_rcv, ysize3, MPI_DOUBLE_PRECISION, south, tag4, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( g_north_snd, ysize3, MPI_DOUBLE_PRECISION, north, tag4, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

! DIAGONALS
 CALL MPI_IRECV(g_ne_rcv, 1, MPI_DOUBLE_PRECISION, ne, tag1, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( g_sw_snd, 1, MPI_DOUBLE_PRECISION, sw, tag1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(g_nw_rcv, 1, MPI_DOUBLE_PRECISION, nw, tag2, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( g_se_snd, 1, MPI_DOUBLE_PRECISION, se, tag2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(g_se_rcv, 1, MPI_DOUBLE_PRECISION, se, tag3, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( g_nw_snd, 1, MPI_DOUBLE_PRECISION, nw, tag3, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(g_sw_rcv, 1, MPI_DOUBLE_PRECISION, sw, tag4, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( g_ne_snd, 1, MPI_DOUBLE_PRECISION, ne, tag4, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)


!--------- Update inward components of g in the real nodes ---------------------
! X DIRECTION
 count = 1
 DO j = yl, yu
   gbar(xu,j,3) = g_east_rcv(count)
   gbar(xu,j,6) = g_east_rcv(count + 1)
   gbar(xu,j,7) = g_east_rcv(count + 2)

   gbar(xl,j,1) = g_west_rcv(count)
   gbar(xl,j,5) = g_west_rcv(count + 1)
   gbar(xl,j,8) = g_west_rcv(count + 2)
   count = count + 3
 END DO

! Y DIRECTION
 count = 1
 DO i = xl, xu
   gbar(i,yu,4) = g_north_rcv(count)
   gbar(i,yu,7) = g_north_rcv(count + 1)
   gbar(i,yu,8) = g_north_rcv(count + 2)

   gbar(i,yl,2) = g_south_rcv(count)
   gbar(i,yl,5) = g_south_rcv(count + 1)
   gbar(i,yl,6) = g_south_rcv(count + 2)
   count = count + 3
 END DO

! DIAGONALS
 gbar(xu,yu,7) = g_ne_rcv
 gbar(xl,yu,8) = g_nw_rcv
 gbar(xl,yl,5) = g_sw_rcv
 gbar(xu,yl,6) = g_se_rcv

 RETURN
 END SUBROUTINE MPI_UpdateG

