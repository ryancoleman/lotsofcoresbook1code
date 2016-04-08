!-------------------------------------------------------------------------------
! Subroutine : PostStreamF
! Revision   : 1.0 (2008-06-15)
! Author     : Carlos Rosales Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
!> @file
!! Update f and phi_f in nodes at processor boundaries using MPI.
!> @details
!! Copy outward f values from the ghosts into inward f values in the real nodes,
!! then update order parameter phi_f and the values values for f in the ghost
!! layers.
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

 SUBROUTINE PostStreamF

! Common Variables
 USE Domain,      ONLY : xl_f, xlg_f, xsize_f, xu_f, xug_f, yl_f, ylg_f, ysize_f, yu_f, yug_f
 USE FluidParams, ONLY : phi_f
 USE LBMParams,   ONLY : f
 USE MPIParams
 USE MPI
 IMPLICIT NONE

! Local Variables
 INTEGER :: count, i, j
 INTEGER :: MPI_ERR, MPI_REQ
 INTEGER :: status(MPI_STATUS_SIZE)

!-------- Save outward pointing values of f ------------------------------------
! X DIRECTION
 count = 1
 DO j = yl_f, yu_f
   f_west_snd(count) = f(xlg_f,j,3)
   f_east_snd(count) = f(xug_f,j,1)
   count = count + 1
 END DO

! Y DIRECTION
 count = 1
 DO i = xl_f, xu_f
   f_south_snd(count) = f(i,ylg_f,4)
   f_north_snd(count) = f(i,yug_f,2)
   count = count + 1
 END DO

!-------- Exchange f values using MPI calls ------------------------------------
! X DIRECTION
 CALL MPI_IRECV(f_east_rcv, xsize_f, MPI_DOUBLE_PRECISION, east, TAG1, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( f_west_snd, xsize_f, MPI_DOUBLE_PRECISION, west, TAG1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(f_west_rcv, xsize_f, MPI_DOUBLE_PRECISION, west, TAG2, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( f_east_snd, xsize_f, MPI_DOUBLE_PRECISION, east, TAG2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

! Y DIRECTION
 CALL MPI_IRECV(f_north_rcv, ysize_f, MPI_DOUBLE_PRECISION, north, TAG1, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( f_south_snd, ysize_f, MPI_DOUBLE_PRECISION, south, TAG1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(f_south_rcv, ysize_f, MPI_DOUBLE_PRECISION, south, TAG2, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( f_north_snd, ysize_f, MPI_DOUBLE_PRECISION, north, TAG2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)


!-------- Update inward components of f in the real nodes ----------------------
! X DIRECTION
 count = 1
 DO j = yl_f, yu_f
   f(xu_f,j,3) = f_east_rcv(count)
   f(xl_f,j,1) = f_west_rcv(count)
   count = count + 1
 END DO

! Y DIRECTION
 count = 1
 DO i = xl_f, xu_f
   f(i,yu_f,4) = f_north_rcv(count)
   f(i,yl_f,2) = f_south_rcv(count)
   count = count + 1
 END DO

!-------- Update the order parameter with the new f values ---------------------
 DO j = yl_f, yu_f
   DO i = xl_f, xu_f
     phi_f(i,j) = f(i,j,0) + f(i,j,1) + f(i,j,2) + f(i,j,3) + f(i,j,4)
   END DO
 END DO

!-------- Exchange and update phi for all ghost nodes --------------------------
! X DIRECTION
 count = 1
 DO j = yl_f, yu_f
   phi_west_snd(count) = phi_f(xl_f,j)
   phi_east_snd(count) = phi_f(xu_f,j)
   count = count + 1
 END DO
 
 CALL MPI_IRECV(phi_east_rcv, xsize_f, MPI_DOUBLE_PRECISION, east, TAG1, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( phi_west_snd, xsize_f, MPI_DOUBLE_PRECISION, west, TAG1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(phi_west_rcv, xsize_f, MPI_DOUBLE_PRECISION, west, TAG2, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( phi_east_snd, xsize_f, MPI_DOUBLE_PRECISION, east, TAG2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)
 
! Y DIRECTION
 count = 1
 DO i = xl_f, xu_f
   phi_south_snd(count) = phi_f(i,yl_f)
   phi_north_snd(count) = phi_f(i,yu_f)
   count = count + 1
 END DO

 CALL MPI_IRECV(phi_north_rcv, ysize_f, MPI_DOUBLE_PRECISION, north, TAG1, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( phi_south_snd, ysize_f, MPI_DOUBLE_PRECISION, south, TAG1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(phi_south_rcv, ysize_f, MPI_DOUBLE_PRECISION, south, TAG2, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( phi_north_snd, ysize_f, MPI_DOUBLE_PRECISION, north, TAG2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

! DIAGONALS
 phi_ne_snd = phi_f(xu_f,yu_f)
 phi_nw_snd = phi_f(xl_f,yu_f)
 phi_sw_snd = phi_f(xl_f,yl_f)
 phi_se_snd = phi_f(xu_f,yl_f)

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

!-------- Update all phi values in the ghost nodes -----------------------------
! X DIRECTION
 count = 1
 DO j = yl_f, yu_f
   phi_f(xug_f,j) = phi_east_rcv(count)
   phi_f(xlg_f,j) = phi_west_rcv(count)
   count = count + 1
 END DO

! Y DIRECTION
 count = 1
 DO i = xl_f, xu_f
   phi_f(i,yug_f) = phi_north_rcv(count)
   phi_f(i,ylg_f) = phi_south_rcv(count)
   count = count + 1
 END DO

! DIAGONALS
 phi_f(xug_f,yug_f) = phi_ne_rcv
 phi_f(xlg_f,yug_f) = phi_nw_rcv
 phi_f(xlg_f,ylg_f) = phi_sw_rcv
 phi_f(xug_f,ylg_f) = phi_se_rcv

 RETURN
 END SUBROUTINE PostStreamF
