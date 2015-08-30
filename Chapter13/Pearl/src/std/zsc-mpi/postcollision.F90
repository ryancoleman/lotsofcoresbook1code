!-------------------------------------------------------------------------------
! Subroutine : PostCollision
! Revision   : 1.0 (2008-06-15)
! Author     : Carlos Rosales Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
!> @file
!! Update g in nodes at processor boundaries and outward f in ghosts using MPI.
!> @details
!! Copy outward g values from the ghosts into inward g values in the real nodes.
!! Copy inward f values from the real nodes into outward f values in the ghosts.
!! These f values are required for the relaxation step in "stream.f".
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

 SUBROUTINE PostCollision

! Common Variables
 USE Domain
 USE LBMParams, ONLY : f, g
 USE MPIParams
 USE MPI
 IMPLICIT NONE

! Local Variables
 INTEGER :: count, i, j
 INTEGER :: MPI_ERR, MPI_REQ
 INTEGER :: status(MPI_STATUS_SIZE)


!-------- Save data for outward f and inward g components ----------------------
! X DIRECTION
 count = 1
 DO j = yl, yu
   f_west_snd(count)          = f(xl ,j,1,now)
   g_west_snd(count)          = g(xlg,j,3,nxt)
   g_west_snd(count + xsize)  = g(xlg,j,6,nxt)
   g_west_snd(count + xsize2) = g(xlg,j,7,nxt)

   f_east_snd(count)          = f(xu ,j,3,now)
   g_east_snd(count)          = g(xug,j,1,nxt)
   g_east_snd(count + xsize)  = g(xug,j,5,nxt)
   g_east_snd(count + xsize2) = g(xug,j,8,nxt)
   count = count + 1
 END DO

! Y DIRECTION
 count = 1
 DO i = xl, xu
   f_south_snd(count)          = f(i,yl ,2,now)
   g_south_snd(count)          = g(i,ylg,4,nxt)
   g_south_snd(count + ysize)  = g(i,ylg,7,nxt)
   g_south_snd(count + ysize2) = g(i,ylg,8,nxt)

   f_north_snd(count)          = f(i,yu ,4,now)
   g_north_snd(count)          = g(i,yug,2,nxt)
   g_north_snd(count + ysize)  = g(i,yug,5,nxt)
   g_north_snd(count + ysize2) = g(i,yug,6,nxt)
   count = count + 1
 END DO

! DIAGONALS
 g_ne_snd = g(xug,yug,5,nxt)
 g_nw_snd = g(xlg,yug,6,nxt)
 g_sw_snd = g(xlg,ylg,7,nxt)
 g_se_snd = g(xug,ylg,8,nxt)

!-------- Exchange information using MPI calls ---------------------------------
! X DIRECTION
 CALL MPI_IRECV(f_east_rcv, xsize, MPI_DOUBLE_PRECISION, east, TAG1, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( f_west_snd, xsize, MPI_DOUBLE_PRECISION, west, TAG1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(g_east_rcv, xsize3, MPI_DOUBLE_PRECISION, east, TAG2, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( g_west_snd, xsize3, MPI_DOUBLE_PRECISION, west, TAG2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(f_west_rcv, xsize, MPI_DOUBLE_PRECISION, west, TAG3, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( f_east_snd, xsize, MPI_DOUBLE_PRECISION, east, TAG3, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(g_west_rcv, xsize3, MPI_DOUBLE_PRECISION, west, TAG4, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( g_east_snd, xsize3, MPI_DOUBLE_PRECISION, east, TAG4, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

! Y DIRECTION
 CALL MPI_IRECV(f_north_rcv, ysize, MPI_DOUBLE_PRECISION, north, TAG1, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( f_south_snd, ysize, MPI_DOUBLE_PRECISION, south, TAG1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(g_north_rcv, ysize3, MPI_DOUBLE_PRECISION, north, TAG2, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( g_south_snd, ysize3, MPI_DOUBLE_PRECISION, south, TAG2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(f_south_rcv, ysize, MPI_DOUBLE_PRECISION, south, TAG3, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( f_north_snd, ysize, MPI_DOUBLE_PRECISION, north, TAG3, MPI_COMM_VGRID, MPI_ERR)
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

!-------- Update outward pointing f components for the ghosts ------------------
! x direction
 count = 1
 DO j = yl, yu
   f(xug,j,1,now) = f_east_rcv(count)
   f(xlg,j,3,now) = f_west_rcv(count)
   count = count + 1
 END DO

! y direction
 count = 1
 DO i = xl, xu
   f(i,yug,2,now) = f_north_rcv(count)
   f(i,ylg,4,now) = f_south_rcv(count)
   count = count + 1
 END DO

 RETURN
 END SUBROUTINE PostCollision

