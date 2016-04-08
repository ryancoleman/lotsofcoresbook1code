!-------------------------------------------------------------------------------
! Subroutine : MPI_UpdateF
! Revision   : 1.0 (2008-06-15)
! Author     : Carlos Rosales Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
!> @file
!! Update fbar in nodes at processor boundaries using MPI.
!> @details
!! Copy outward values of the order parameter distribution function fbar from
!! the ghosts into inward values in the real nodes after the streaming step.
!!
!! MPI exchange is done using non-blocking receive and blocking send pairs
!! (completed by a wait statement) in alternate directions to avoid deadlock.
!! The components of fbar that need to be exchanged are packed into buffer
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
!--------------------------------------------------------------------------------

 SUBROUTINE MPI_UpdateF

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

!--------- Save data for inward f components -----------------------------------
! X DIRECTION
 count = 1
 DO j = yl_f, yu_f
   f_west_snd(count)     = fbar(xlg_f,j,3)
   f_west_snd(count + 1) = fbar(xlg_f,j,6)
   f_west_snd(count + 2) = fbar(xlg_f,j,7)

   f_east_snd(count)     = fbar(xug_f,j,1)
   f_east_snd(count + 1) = fbar(xug_f,j,5)
   f_east_snd(count + 2) = fbar(xug_f,j,8)
   count = count + 3
 END DO

! Y DIRECTION
 count = 1
 DO i = xl_f, xu_f
   f_south_snd(count)     = fbar(i,ylg_f,4)
   f_south_snd(count + 1) = fbar(i,ylg_f,7)
   f_south_snd(count + 2) = fbar(i,ylg_f,8)

   f_north_snd(count)     = fbar(i,yug_f,2)
   f_north_snd(count + 1) = fbar(i,yug_f,5)
   f_north_snd(count + 2) = fbar(i,yug_f,6)
   count = count + 3
 END DO

! DIAGONALS
 f_ne_snd = fbar(xug_f,yug_f,5)
 f_nw_snd = fbar(xlg_f,yug_f,6)
 f_sw_snd = fbar(xlg_f,ylg_f,7)
 f_se_snd = fbar(xug_f,ylg_f,8)

!--------- Exchange information using MPI calls --------------------------------
! X DIRECTION
 CALL MPI_IRECV(f_east_rcv, xsize3_f, MPI_DOUBLE_PRECISION, east, tag1, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( f_west_snd, xsize3_f, MPI_DOUBLE_PRECISION, west, tag1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(f_west_rcv, xsize3_f, MPI_DOUBLE_PRECISION, west, tag3, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( f_east_snd, xsize3_f, MPI_DOUBLE_PRECISION, east, tag3, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

! Y DIRECTION
 CALL MPI_IRECV(f_north_rcv, ysize3_f, MPI_DOUBLE_PRECISION, north, tag1, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( f_south_snd, ysize3_f, MPI_DOUBLE_PRECISION, south, tag1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(f_south_rcv, ysize3_f, MPI_DOUBLE_PRECISION, south, tag3, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( f_north_snd, ysize3_f, MPI_DOUBLE_PRECISION, north, tag3, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

! DIAGONALS
 CALL MPI_IRECV(f_ne_rcv, 1, MPI_DOUBLE_PRECISION, ne, tag1, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( f_sw_snd, 1, MPI_DOUBLE_PRECISION, sw, tag1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(f_nw_rcv, 1, MPI_DOUBLE_PRECISION, nw, tag2, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( f_se_snd, 1, MPI_DOUBLE_PRECISION, se, tag2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(f_se_rcv, 1, MPI_DOUBLE_PRECISION, se, tag3, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( f_nw_snd, 1, MPI_DOUBLE_PRECISION, nw, tag3, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(f_sw_rcv, 1, MPI_DOUBLE_PRECISION, sw, tag4, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( f_ne_snd, 1, MPI_DOUBLE_PRECISION, ne, tag4, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)


!--------- Update inward components of f in the real nodes ---------------------
! X DIRECTION
 count = 1
 DO j = yl_f, yu_f
   fbar(xu_f,j,3) = f_east_rcv(count)
   fbar(xu_f,j,6) = f_east_rcv(count + 1)
   fbar(xu_f,j,7) = f_east_rcv(count + 2)

   fbar(xl_f,j,1) = f_west_rcv(count)
   fbar(xl_f,j,5) = f_west_rcv(count + 1)
   fbar(xl_f,j,8) = f_west_rcv(count + 2)
   count = count + 3
 END DO

! Y DIRECTION
 count = 1
 DO i = xl_f, xu_f
   fbar(i,yu_f,4) = f_north_rcv(count)
   fbar(i,yu_f,7) = f_north_rcv(count + 1) 
   fbar(i,yu_f,8) = f_north_rcv(count + 2)

   fbar(i,yl_f,2) = f_south_rcv(count)
   fbar(i,yl_f,5) = f_south_rcv(count + 1)
   fbar(i,yl_f,6) = f_south_rcv(count + 2)
   count = count + 3
 END DO

! DIAGONALS
 fbar(xu_f,yu_f,7) = f_ne_rcv
 fbar(xl_f,yu_f,8) = f_nw_rcv
 fbar(xl_f,yl_f,5) = f_sw_rcv
 fbar(xu_f,yl_f,6) = f_se_rcv


 RETURN
 END SUBROUTINE MPI_UpdateF

