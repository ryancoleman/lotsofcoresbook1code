!-------------------------------------------------------------------------------
! Subroutine : MPI_UpdateRho
! Revision   : 1.0 (2008-06-15)
! Author     : Carlos Rosales Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
!> @file
!! Update rho values in the three ghost layers using MPI.
!> @details
!! Copy outward values of the density (rho) from the three outermost real node
!! layers onto the three ghost layers.
!!
!! MPI exchange is done using non-blocking receive and blocking send pairs
!! (completed by a wait statement) in alternate directions to avoid deadlock.
!! Rho is packed into arrays for each of the directions before the MPI exchange.

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

 SUBROUTINE MPI_UpdateRho

! Common Variables
 USE Domain
 USE FluidParams, ONLY : rho
 USE MPIParams
 USE MPI
 IMPLICIT NONE

! Local Variables
 INTEGER :: count, i, j
 INTEGER :: MPI_ERR, MPI_REQ
 INTEGER :: status(MPI_STATUS_SIZE)

!--------- Save data from real nodes -------------------------------------------
! X DIRECTION
 count = 1
 DO j = yl, yu
   rho_west_snd(count)     = rho(xl  ,j)
   rho_west_snd(count + 1) = rho(xl+1,j)
   rho_west_snd(count + 2) = rho(xl+2,j)
   rho_east_snd(count)     = rho(xu  ,j)
   rho_east_snd(count + 1) = rho(xu-1,j)
   rho_east_snd(count + 2) = rho(xu-2,j)
   count = count + 3
 END DO

! Y DIRECTION
 count = 1
 DO i = xl, xu
   rho_south_snd(count)     = rho(i,yl  )
   rho_south_snd(count + 1) = rho(i,yl+1)
   rho_south_snd(count + 2) = rho(i,yl+2)
   rho_north_snd(count)     = rho(i,yu  )
   rho_north_snd(count + 1) = rho(i,yu-1)
   rho_north_snd(count + 2) = rho(i,yu-2)
   count = count + 3
 END DO

! NE DIAGONAL
 rho_ne_snd(1) = rho(xu  ,yu  )
 rho_ne_snd(2) = rho(xu-1,yu  )
 rho_ne_snd(3) = rho(xu-2,yu  )
 rho_ne_snd(4) = rho(xu  ,yu-1)
 rho_ne_snd(5) = rho(xu-1,yu-1)
 rho_ne_snd(6) = rho(xu-2,yu-1)
 rho_ne_snd(7) = rho(xu  ,yu-2)
 rho_ne_snd(8) = rho(xu-1,yu-2)
 rho_ne_snd(9) = rho(xu-2,yu-2)
 
! NW DIAGONAL
 rho_nw_snd(1) = rho(xl  ,yu  )
 rho_nw_snd(2) = rho(xl+1,yu  )
 rho_nw_snd(3) = rho(xl+2,yu  )
 rho_nw_snd(4) = rho(xl  ,yu-1)
 rho_nw_snd(5) = rho(xl+1,yu-1)
 rho_nw_snd(6) = rho(xl+2,yu-1)
 rho_nw_snd(7) = rho(xl  ,yu-2)
 rho_nw_snd(8) = rho(xl+1,yu-2)
 rho_nw_snd(9) = rho(xl+2,yu-2)

! SW DIAGONAL
 rho_sw_snd(1) = rho(xl  ,yl  )
 rho_sw_snd(2) = rho(xl+1,yl  )
 rho_sw_snd(3) = rho(xl+2,yl  )
 rho_sw_snd(4) = rho(xl  ,yl+1)
 rho_sw_snd(5) = rho(xl+1,yl+1)
 rho_sw_snd(6) = rho(xl+2,yl+1)
 rho_sw_snd(7) = rho(xl  ,yl+2)
 rho_sw_snd(8) = rho(xl+1,yl+2)
 rho_sw_snd(9) = rho(xl+2,yl+2)

! SE DIAGONAL 
 rho_se_snd(1) = rho(xu  ,yl  ) 
 rho_se_snd(2) = rho(xu-1,yl  )
 rho_se_snd(3) = rho(xu-2,yl  )
 rho_se_snd(4) = rho(xu  ,yl+1) 
 rho_se_snd(5) = rho(xu-1,yl+1)
 rho_se_snd(6) = rho(xu-2,yl+1)
 rho_se_snd(7) = rho(xu  ,yl+2) 
 rho_se_snd(8) = rho(xu-1,yl+2)
 rho_se_snd(9) = rho(xu-2,yl+2)

 
!--------- Exchange information using MPI calls --------------------------------
! X DIRECTION
 CALL MPI_IRECV(rho_east_rcv, xsize3, MPI_DOUBLE_PRECISION, east, TAG1, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( rho_west_snd, xsize3, MPI_DOUBLE_PRECISION, west, TAG1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(rho_west_rcv, xsize3, MPI_DOUBLE_PRECISION, west, TAG2, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( rho_east_snd, xsize3, MPI_DOUBLE_PRECISION, east, TAG2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

! Y DIRECTION
 CALL MPI_IRECV(rho_north_rcv, ysize3, MPI_DOUBLE_PRECISION, north, TAG3, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( rho_south_snd, ysize3, MPI_DOUBLE_PRECISION, south, TAG3, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(rho_south_rcv, ysize3, MPI_DOUBLE_PRECISION, south, TAG4, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( rho_north_snd, ysize3, MPI_DOUBLE_PRECISION, north, TAG4, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

! DIAGONALS
 CALL MPI_IRECV(rho_ne_rcv, 9, MPI_DOUBLE_PRECISION, ne, TAG1, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( rho_sw_snd, 9, MPI_DOUBLE_PRECISION, sw, TAG1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(rho_nw_rcv, 9, MPI_DOUBLE_PRECISION, nw, TAG2, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( rho_se_snd, 9, MPI_DOUBLE_PRECISION, se, TAG2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(rho_se_rcv, 9, MPI_DOUBLE_PRECISION, se, TAG3, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( rho_nw_snd, 9, MPI_DOUBLE_PRECISION, nw, TAG3, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(rho_sw_rcv, 9, MPI_DOUBLE_PRECISION, sw, TAG4, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( rho_ne_snd, 9, MPI_DOUBLE_PRECISION, ne, TAG4, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)
 

!---------  Update the density rho in the three ghost layers -------------------
! X DIRECTION
 count = 1
 DO j = yl, yu
   rho(xug ,j) = rho_east_rcv(count)
   rho(xug2,j) = rho_east_rcv(count + 1)
   rho(xug3,j) = rho_east_rcv(count + 2)
   rho(xlg ,j) = rho_west_rcv(count)
   rho(xlg2,j) = rho_west_rcv(count + 1)
   rho(xlg3,j) = rho_west_rcv(count + 2)
   count = count + 3
 END DO

! Y DIRECTION
 count = 1
 DO i = xl, xu
   rho(i,yug ) = rho_north_rcv(count)
   rho(i,yug2) = rho_north_rcv(count + 1)
   rho(i,yug3) = rho_north_rcv(count + 2)
   rho(i,ylg ) = rho_south_rcv(count)
   rho(i,ylg2) = rho_south_rcv(count + 1)
   rho(i,ylg3) = rho_south_rcv(count + 2)
   count = count + 3
 END DO

! NE DIAGONAL
 rho(xug ,yug ) = rho_ne_rcv(1)
 rho(xug2,yug ) = rho_ne_rcv(2)
 rho(xug3,yug ) = rho_ne_rcv(3)
 rho(xug ,yug2) = rho_ne_rcv(4)
 rho(xug2,yug2) = rho_ne_rcv(5)
 rho(xug3,yug2) = rho_ne_rcv(6)
 rho(xug ,yug3) = rho_ne_rcv(7)
 rho(xug2,yug3) = rho_ne_rcv(8)
 rho(xug3,yug3) = rho_ne_rcv(9)

! NW DIAGONAL 
 rho(xlg ,yug ) = rho_nw_rcv(1)
 rho(xlg2,yug ) = rho_nw_rcv(2)
 rho(xlg3,yug ) = rho_nw_rcv(3)
 rho(xlg ,yug2) = rho_nw_rcv(4)
 rho(xlg2,yug2) = rho_nw_rcv(5)
 rho(xlg3,yug2) = rho_nw_rcv(6)
 rho(xlg ,yug3) = rho_nw_rcv(7)
 rho(xlg2,yug3) = rho_nw_rcv(8)
 rho(xlg3,yug3) = rho_nw_rcv(9)

! SW DIAGONAL 
 rho(xlg ,ylg ) = rho_sw_rcv(1)
 rho(xlg2,ylg ) = rho_sw_rcv(2)
 rho(xlg3,ylg ) = rho_sw_rcv(3)
 rho(xlg ,ylg2) = rho_sw_rcv(4)
 rho(xlg2,ylg2) = rho_sw_rcv(5)
 rho(xlg3,ylg2) = rho_sw_rcv(6)
 rho(xlg ,ylg3) = rho_sw_rcv(7)
 rho(xlg2,ylg3) = rho_sw_rcv(8)
 rho(xlg3,ylg3) = rho_sw_rcv(9)

! SE DIAGONAL
 rho(xug ,ylg ) = rho_se_rcv(1)
 rho(xug2,ylg ) = rho_se_rcv(2)
 rho(xug3,ylg ) = rho_se_rcv(3)
 rho(xug ,ylg2) = rho_se_rcv(4)
 rho(xug2,ylg2) = rho_se_rcv(5)
 rho(xug3,ylg2) = rho_se_rcv(6)
 rho(xug ,ylg3) = rho_se_rcv(7)
 rho(xug2,ylg3) = rho_se_rcv(8)
 rho(xug3,ylg3) = rho_se_rcv(9)

 RETURN
 END SUBROUTINE MPI_UpdateRho

