!-------------------------------------------------------------------------------
! Subroutine : MPI_UpdateRhoF
! Revision   : 1.0 (2008-06-15)
! Author     : Carlos Rosales Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
!> @file
!! Update rho_f values in the three ghost layers using MPI.
!> @details
!! Copy outward values of the density (rho_f) from the three outermost real node
!! layers onto the three ghost layers for the order parameter grid.
!!
!! MPI exchange is done using non-blocking receive and blocking send pairs
!! (completed by a wait statement) in alternate directions to avoid deadlock.
!! Rho_f is packed into arrays for each of the directions before the MPI
!! exchange.

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

 SUBROUTINE MPI_UpdateRhoF

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


!--------- Calculate the updated density rho_f ---------------------------------
 DO j = yl_f, yu_f
   DO i = xl_f, xu_f
      rho_f(i,j) = fbar(i,j,0) + fbar(i,j,1) + fbar(i,j,2) + fbar(i,j,3) + fbar(i,j,4) &
                 + fbar(i,j,5) + fbar(i,j,6) + fbar(i,j,7) + fbar(i,j,8)
   END DO
 END DO

!--------- Save data from real nodes -------------------------------------------
! X DIRECTION
 count = 1
 DO j = yl_f, yu_f
   rho_west_snd(count)     = rho_f(xl_f  ,j)
   rho_west_snd(count + 1) = rho_f(xl_f+1,j)
   rho_west_snd(count + 2) = rho_f(xl_f+2,j)

   rho_east_snd(count)     = rho_f(xu_f  ,j)
   rho_east_snd(count + 1) = rho_f(xu_f-1,j)
   rho_east_snd(count + 2) = rho_f(xu_f-2,j)
   count = count + 3
 END DO

! Y DIRECTION
 count = 1
 DO i = xl_f, xu_f
   rho_south_snd(count)     = rho_f(i,yl_f  )
   rho_south_snd(count + 1) = rho_f(i,yl_f+1)
   rho_south_snd(count + 2) = rho_f(i,yl_f+2)

   rho_north_snd(count)     = rho_f(i,yu_f  )
   rho_north_snd(count + 1) = rho_f(i,yu_f-1)
   rho_north_snd(count + 2) = rho_f(i,yu_f-2)
   count = count + 3
 END DO

! NE DIAGONAL
 rho_ne_snd(1) = rho_f(xu_f  ,yu_f  )
 rho_ne_snd(2) = rho_f(xu_f-1,yu_f  )
 rho_ne_snd(3) = rho_f(xu_f-2,yu_f  )
 rho_ne_snd(4) = rho_f(xu_f  ,yu_f-1)
 rho_ne_snd(5) = rho_f(xu_f-1,yu_f-1)
 rho_ne_snd(6) = rho_f(xu_f-2,yu_f-1)
 rho_ne_snd(7) = rho_f(xu_f  ,yu_f-2)
 rho_ne_snd(8) = rho_f(xu_f-1,yu_f-2)
 rho_ne_snd(9) = rho_f(xu_f-2,yu_f-2)
 
! NW DIAGONAL
 rho_nw_snd(1) = rho_f(xl_f  ,yu_f  )
 rho_nw_snd(2) = rho_f(xl_f+1,yu_f  )
 rho_nw_snd(3) = rho_f(xl_f+2,yu_f  )
 rho_nw_snd(4) = rho_f(xl_f  ,yu_f-1)
 rho_nw_snd(5) = rho_f(xl_f+1,yu_f-1)
 rho_nw_snd(6) = rho_f(xl_f+2,yu_f-1)
 rho_nw_snd(7) = rho_f(xl_f  ,yu_f-2)
 rho_nw_snd(8) = rho_f(xl_f+1,yu_f-2)
 rho_nw_snd(9) = rho_f(xl_f+2,yu_f-2)

! SW DIAGONAL
 rho_sw_snd(1) = rho_f(xl_f  ,yl_f  )
 rho_sw_snd(2) = rho_f(xl_f+1,yl_f  )
 rho_sw_snd(3) = rho_f(xl_f+2,yl_f  )
 rho_sw_snd(4) = rho_f(xl_f  ,yl_f+1)
 rho_sw_snd(5) = rho_f(xl_f+1,yl_f+1)
 rho_sw_snd(6) = rho_f(xl_f+2,yl_f+1)
 rho_sw_snd(7) = rho_f(xl_f  ,yl_f+2)
 rho_sw_snd(8) = rho_f(xl_f+1,yl_f+2)
 rho_sw_snd(9) = rho_f(xl_f+2,yl_f+2)

! SE DIAGONAL 
 rho_se_snd(1) = rho_f(xu_f  ,yl_f  )
 rho_se_snd(2) = rho_f(xu_f-1,yl_f  )
 rho_se_snd(3) = rho_f(xu_f-2,yl_f  )
 rho_se_snd(4) = rho_f(xu_f  ,yl_f+1)
 rho_se_snd(5) = rho_f(xu_f-1,yl_f+1)
 rho_se_snd(6) = rho_f(xu_f-2,yl_f+1)
 rho_se_snd(7) = rho_f(xu_f  ,yl_f+2)
 rho_se_snd(8) = rho_f(xu_f-1,yl_f+2)
 rho_se_snd(9) = rho_f(xu_f-2,yl_f+2)


!--------- Exchange information using MPI calls --------------------------------
! X DIRECTION
 CALL MPI_IRECV(rho_east_rcv, xsize3_f, MPI_DOUBLE_PRECISION, east, tag1, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( rho_west_snd, xsize3_f, MPI_DOUBLE_PRECISION, west, tag1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(rho_west_rcv, xsize3_f, MPI_DOUBLE_PRECISION, west, tag2, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( rho_east_snd, xsize3_f, MPI_DOUBLE_PRECISION, east, tag2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

! Y DIRECTION
 CALL MPI_IRECV(rho_north_rcv, ysize3_f, MPI_DOUBLE_PRECISION, north, tag3, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( rho_south_snd, ysize3_f, MPI_DOUBLE_PRECISION, south, tag3, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(rho_south_rcv, ysize3_f, MPI_DOUBLE_PRECISION, south, tag4, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( rho_north_snd, ysize3_f, MPI_DOUBLE_PRECISION, north, tag4, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

! DIAGONALS
 CALL MPI_IRECV(rho_ne_rcv, 9, MPI_DOUBLE_PRECISION, ne, tag1, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( rho_sw_snd, 9, MPI_DOUBLE_PRECISION, sw, tag1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(rho_nw_rcv, 9, MPI_DOUBLE_PRECISION, nw, tag2, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( rho_se_snd, 9, MPI_DOUBLE_PRECISION, se, tag2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(rho_se_rcv, 9, MPI_DOUBLE_PRECISION, se, tag3, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( rho_nw_snd, 9, MPI_DOUBLE_PRECISION, nw, tag3, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(rho_sw_rcv, 9, MPI_DOUBLE_PRECISION, sw, tag4, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( rho_ne_snd, 9, MPI_DOUBLE_PRECISION, ne, tag4, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)


!--------- Update the density rho_f in the three ghost layers  -----------------
! X DIRECTION
 count = 1
 DO j = yl_f, yu_f
   rho_f(xug_f ,j) = rho_east_rcv(count)
   rho_f(xug2_f,j) = rho_east_rcv(count + 1)
   rho_f(xug3_f,j) = rho_east_rcv(count + 2)

   rho_f(xlg_f ,j) = rho_west_rcv(count)
   rho_f(xlg2_f,j) = rho_west_rcv(count + 1)
   rho_f(xlg3_f,j) = rho_west_rcv(count + 2)
   count = count + 3
 END DO

! Y DIRECTION
 count = 1
 DO i = xl_f, xu_f
   rho_f(i,yug_f ) = rho_north_rcv(count)
   rho_f(i,yug2_f) = rho_north_rcv(count + 1)
   rho_f(i,yug3_f) = rho_north_rcv(count + 2)

   rho_f(i,ylg_f ) = rho_south_rcv(count)
   rho_f(i,ylg2_f) = rho_south_rcv(count + 1)
   rho_f(i,ylg3_f) = rho_south_rcv(count + 2)
   count = count + 3
 END DO

! NE DIAGONAL
 rho_f(xug_f ,yug_f ) = rho_ne_rcv(1)
 rho_f(xug2_f,yug_f ) = rho_ne_rcv(2)
 rho_f(xug3_f,yug_f ) = rho_ne_rcv(3)
 rho_f(xug_f ,yug2_f) = rho_ne_rcv(4)
 rho_f(xug2_f,yug2_f) = rho_ne_rcv(5)
 rho_f(xug3_f,yug2_f) = rho_ne_rcv(6)
 rho_f(xug_f ,yug3_f) = rho_ne_rcv(7)
 rho_f(xug2_f,yug3_f) = rho_ne_rcv(8)
 rho_f(xug3_f,yug3_f) = rho_ne_rcv(9)

! NW DIAGONAL 
 rho_f(xlg_f ,yug_f ) = rho_nw_rcv(1)
 rho_f(xlg2_f,yug_f ) = rho_nw_rcv(2)
 rho_f(xlg3_f,yug_f ) = rho_nw_rcv(3)
 rho_f(xlg_f ,yug2_f) = rho_nw_rcv(4)
 rho_f(xlg2_f,yug2_f) = rho_nw_rcv(5)
 rho_f(xlg3_f,yug2_f) = rho_nw_rcv(6)
 rho_f(xlg_f ,yug3_f) = rho_nw_rcv(7)
 rho_f(xlg2_f,yug3_f) = rho_nw_rcv(8)
 rho_f(xlg3_f,yug3_f) = rho_nw_rcv(9)

! SW DIAGONAL 
 rho_f(xlg_f ,ylg_f ) = rho_sw_rcv(1)
 rho_f(xlg2_f,ylg_f ) = rho_sw_rcv(2)
 rho_f(xlg3_f,ylg_f ) = rho_sw_rcv(3)
 rho_f(xlg_f ,ylg2_f) = rho_sw_rcv(4)
 rho_f(xlg2_f,ylg2_f) = rho_sw_rcv(5)
 rho_f(xlg3_f,ylg2_f) = rho_sw_rcv(6)
 rho_f(xlg_f ,ylg3_f) = rho_sw_rcv(7)
 rho_f(xlg2_f,ylg3_f) = rho_sw_rcv(8)
 rho_f(xlg3_f,ylg3_f) = rho_sw_rcv(9)

! SE DIAGONAL
 rho_f(xug_f ,ylg_f ) = rho_se_rcv(1)
 rho_f(xug2_f,ylg_f ) = rho_se_rcv(2)
 rho_f(xug3_f,ylg_f ) = rho_se_rcv(3)
 rho_f(xug_f ,ylg2_f) = rho_se_rcv(4)
 rho_f(xug2_f,ylg2_f) = rho_se_rcv(5)
 rho_f(xug3_f,ylg2_f) = rho_se_rcv(6)
 rho_f(xug_f ,ylg3_f) = rho_se_rcv(7)
 rho_f(xug2_f,ylg3_f) = rho_se_rcv(8)
 rho_f(xug3_f,ylg3_f) = rho_se_rcv(9)

 RETURN
 END SUBROUTINE MPI_UpdateRhoF

