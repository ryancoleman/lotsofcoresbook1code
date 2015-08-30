!-------------------------------------------------------------------------------
! Subroutine : PostStream
! Revision   : 1.0 (2008-06-15)
! Author     : Carlos Rosales Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
!> @file
!! Update f in nodes at processor boundaries using MPI.
!> @details
!! Copy outward f values from the ghosts into inward f values in the real nodes,
!! then update order parameter phi and the values values for g and f in the
!! ghost layer.
!!
!! MPI exchange is done using non-blocking receive and blocking send pairs
!! (completed by a wait statement) in alternate directions to avoid deadlock.
!! The componets of g and f that need to be exchanged are packed into buffer
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

 SUBROUTINE PostStream

! Common Variables
 USE Domain
 USE FluidParams, ONLY : phi
 USE LBMParams
 USE MPIParams
 USE MPI
 IMPLICIT NONE

! Local Variables
 INTEGER :: count, i, j
 INTEGER :: MPI_ERR, MPI_REQ
 INTEGER :: status(MPI_STATUS_SIZE)

!-------- Exchange inward pointing values of f ---------------------------------
! X DIRECTION
 count = 1
 DO j = yl, yu
   f_west_snd(count) = f(xlg,j,3,nxt)
   f_east_snd(count) = f(xug,j,1,nxt)
   count = count + 1
 END DO

 CALL MPI_IRECV(f_east_rcv, xsize, MPI_DOUBLE_PRECISION, east, TAG1, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( f_west_snd, xsize, MPI_DOUBLE_PRECISION, west, TAG1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(f_west_rcv, xsize, MPI_DOUBLE_PRECISION, west, TAG2, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( f_east_snd, xsize, MPI_DOUBLE_PRECISION, east, TAG2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

! Y DIRECTION
 count = 1
 DO i = xl, xu
   f_south_snd(count) = f(i,ylg,4,nxt)
   f_north_snd(count) = f(i,yug,2,nxt)
   count = count + 1
 END DO
 
 CALL MPI_IRECV(f_north_rcv, ysize, MPI_DOUBLE_PRECISION, north, TAG1, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( f_south_snd, ysize, MPI_DOUBLE_PRECISION, south, TAG1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(f_south_rcv, ysize, MPI_DOUBLE_PRECISION, south, TAG2, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( f_north_snd, ysize, MPI_DOUBLE_PRECISION, north, TAG2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)


!-------- Update inward components of f and g in the real nodes ----------------
! X DIRECTION
 count = 1
 DO j = yl, yu
   f(xu,j,3,nxt) = f_east_rcv(count)
   g(xu,j,3,nxt) = g_east_rcv(count)
   g(xu,j,6,nxt) = g_east_rcv(count + xsize)
   g(xu,j,7,nxt) = g_east_rcv(count + xsize2)

   f(xl,j,1,nxt) = f_west_rcv(count)
   g(xl,j,1,nxt) = g_west_rcv(count)
   g(xl,j,5,nxt) = g_west_rcv(count + xsize)
   g(xl,j,8,nxt) = g_west_rcv(count + xsize2)
   count = count + 1
 END DO

! Y DIRECTION
 count = 1
 DO i = xl, xu
   f(i,yu,4,nxt) = f_north_rcv(count)
   g(i,yu,4,nxt) = g_north_rcv(count)
   g(i,yu,7,nxt) = g_north_rcv(count + ysize)
   g(i,yu,8,nxt) = g_north_rcv(count + ysize2)

   f(i,yl,2,nxt) = f_south_rcv(count)
   g(i,yl,2,nxt) = g_south_rcv(count)
   g(i,yl,5,nxt) = g_south_rcv(count + ysize)
   g(i,yl,6,nxt) = g_south_rcv(count + ysize2)
   count = count + 1
 END DO

! DIAGONALS
 g(xu,yu,7,nxt) = g_ne_rcv
 g(xl,yu,8,nxt) = g_nw_rcv
 g(xl,yl,5,nxt) = g_sw_rcv
 g(xu,yl,6,nxt) = g_se_rcv

!-------- Calculate order parameter using updated values of f ------------------
 DO j = yl, yu
   DO i = xl, xu
     phi(i,j) = f(i,j,0,nxt) + f(i,j,1,nxt) + f(i,j,2,nxt) + f(i,j,3,nxt) &
              + f(i,j,4,nxt)
   END DO
 END DO

!-------- Exchange and update phi for all ghost nodes --------------------------
! X DIRECTION
 count = 1
 DO j = yl, yu
   phi_west_snd(count) = phi(xl,j)
   phi_east_snd(count) = phi(xu,j)
   count = count + 1
 END DO

 CALL MPI_IRECV(phi_east_rcv, xsize, MPI_DOUBLE_PRECISION, east, TAG1, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( phi_west_snd, xsize, MPI_DOUBLE_PRECISION, west, TAG1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(phi_west_rcv, xsize, MPI_DOUBLE_PRECISION, west, TAG2, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( phi_east_snd, xsize, MPI_DOUBLE_PRECISION, east, TAG2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

! Y DIRECTION
 count = 1
 DO i = xl, xu
   phi_south_snd(count) = phi(i,yl)
   phi_north_snd(count) = phi(i,yu)
   count = count + 1
 END DO

 CALL MPI_IRECV(phi_north_rcv, ysize, MPI_DOUBLE_PRECISION, north, TAG1, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( phi_south_snd, ysize, MPI_DOUBLE_PRECISION, south, TAG1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(phi_south_rcv, ysize, MPI_DOUBLE_PRECISION, south, TAG2, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( phi_north_snd, ysize, MPI_DOUBLE_PRECISION, north, TAG2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

! DIAGONALS
 phi_ne_snd = phi(xu,yu)
 phi_nw_snd = phi(xl,yu)
 phi_sw_snd = phi(xl,yl)
 phi_se_snd = phi(xu,yl)

 CALL MPI_IRECV(phi_ne_rcv, 1, MPI_DOUBLE_PRECISION, ne, TAG1, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( phi_sw_snd, 1, MPI_DOUBLE_PRECISION, sw, TAG1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(phi_sw_rcv, 1, MPI_DOUBLE_PRECISION, sw, TAG2, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( phi_ne_snd, 1, MPI_DOUBLE_PRECISION, ne, TAG2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(phi_se_rcv, 1, MPI_DOUBLE_PRECISION, se, TAG1, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( phi_nw_snd, 1, MPI_DOUBLE_PRECISION, nw, TAG1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(phi_nw_rcv, 1, MPI_DOUBLE_PRECISION, nw, TAG2, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( phi_se_snd, 1, MPI_DOUBLE_PRECISION, se, TAG2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

!---------Update all phi values in the ghost nodes---------------------------------------------------------------------
! X DIRECTION
 count = 1
 DO j = yl, yu
   phi(xug,j) = phi_east_rcv(count)
   phi(xlg,j) = phi_west_rcv(count)
   count = count + 1
 END DO

! Y DIRECTION
 count = 1
 DO i = xl, xu
   phi(i,yug) = phi_north_rcv(count)
   phi(i,ylg) = phi_south_rcv(count)
   count = count + 1
 END DO

! DIAGONALS
 phi(xug,yug) = phi_ne_rcv
 phi(xlg,yug) = phi_nw_rcv
 phi(xlg,ylg) = phi_sw_rcv
 phi(xug,ylg) = phi_se_rcv

 RETURN
 END SUBROUTINE PostStream
