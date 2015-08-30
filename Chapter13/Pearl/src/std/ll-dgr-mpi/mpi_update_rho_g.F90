!-------------------------------------------------------------------------------
! Subroutine : MPI_UpdateRhoG
! Revision   : 1.0 (2008-06-15)
! Author     : Carlos Rosales Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
!> @file
!! Update rho values in the two ghost layers using MPI.
!> @details
!! Copy outward values of the density (rho) from the two outermost real node
!! layers onto the two ghost layers for the pressure distribution grid.
!!
!! MPI exchange is done using non-blocking receive and blocking send pairs
!! (completed by a wait statement) in alternate directions to avoid deadlock.
!! Rho is packed into arrays for each of the directions before the MPI
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

 SUBROUTINE MPI_UpdateRhoG

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


!--------- Save data from real nodes -------------------------------------------
! X DIRECTION
 count = 1
 DO j = yl, yu
   rhoG_west_snd(count)     = rho(xl  ,j)
   rhoG_west_snd(count + 1) = rho(xl+1,j)
   rhoG_east_snd(count)     = rho(xu  ,j)
   rhoG_east_snd(count + 1) = rho(xu-1,j)
   count = count + 2
 END DO
 
! Y DIRECTION
 count = 1
 DO i = xl, xu
   rhoG_south_snd(count)     = rho(i,yl  )
   rhoG_south_snd(count + 1) = rho(i,yl+1)
   rhoG_north_snd(count)     = rho(i,yu  )
   rhoG_north_snd(count + 1) = rho(i,yu-1)
   count = count + 2
 END DO
 
! NE DIAGONAL
 rhoG_ne_snd(1) = rho(xu  ,yu  )
 rhoG_ne_snd(2) = rho(xu-1,yu  )
 rhoG_ne_snd(3) = rho(xu  ,yu-1)
 rhoG_ne_snd(4) = rho(xu-1,yu-1)

! NW DIAGONAL
 rhoG_nw_snd(1) = rho(xl  ,yu  )
 rhoG_nw_snd(2) = rho(xl+1,yu  )
 rhoG_nw_snd(3) = rho(xl  ,yu-1)
 rhoG_nw_snd(4) = rho(xl+1,yu-1)

! SW DIAGONAL
 rhoG_sw_snd(1) = rho(xl  ,yl  )
 rhoG_sw_snd(2) = rho(xl+1,yl  )
 rhoG_sw_snd(3) = rho(xl  ,yl+1)
 rhoG_sw_snd(4) = rho(xl+1,yl+1)

! SE DIAGONAL 
 rhoG_se_snd(1) = rho(xu  ,yl  )
 rhoG_se_snd(2) = rho(xu-1,yl  )
 rhoG_se_snd(3) = rho(xu  ,yl+1)
 rhoG_se_snd(4) = rho(xu-1,yl+1)


!--------- Exchange information using MPI calls --------------------------------
! X DIRECTION
 CALL MPI_IRECV(rhoG_east_rcv, xsize2, MPI_DOUBLE_PRECISION, east, tag1, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( rhoG_west_snd, xsize2, MPI_DOUBLE_PRECISION, west, tag1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(rhoG_west_rcv, xsize2, MPI_DOUBLE_PRECISION, west, tag2, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( rhoG_east_snd, xsize2, MPI_DOUBLE_PRECISION, east, tag2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

! Y DIRECTION
 CALL MPI_IRECV(rhoG_north_rcv, ysize2, MPI_DOUBLE_PRECISION, north, tag3, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( rhoG_south_snd, ysize2, MPI_DOUBLE_PRECISION, south, tag3, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(rhoG_south_rcv, ysize2, MPI_DOUBLE_PRECISION, south, tag4, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( rhoG_north_snd, ysize2, MPI_DOUBLE_PRECISION, north, tag4, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

! DIAGONALS
 CALL MPI_IRECV(rhoG_ne_rcv, 4, MPI_DOUBLE_PRECISION, ne, tag1, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( rhoG_sw_snd, 4, MPI_DOUBLE_PRECISION, sw, tag1, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(rhoG_nw_rcv, 4, MPI_DOUBLE_PRECISION, nw, tag2, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( rhoG_se_snd, 4, MPI_DOUBLE_PRECISION, se, tag2, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(rhoG_se_rcv, 4, MPI_DOUBLE_PRECISION, se, tag3, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( rhoG_nw_snd, 4, MPI_DOUBLE_PRECISION, nw, tag3, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)

 CALL MPI_IRECV(rhoG_sw_rcv, 4, MPI_DOUBLE_PRECISION, sw, tag4, MPI_COMM_VGRID, MPI_REQ, MPI_ERR)
 CALL MPI_SEND( rhoG_ne_snd, 4, MPI_DOUBLE_PRECISION, ne, tag4, MPI_COMM_VGRID, MPI_ERR)
 CALL MPI_WAIT(MPI_REQ, status, MPI_ERR)


!--------- Update the density rho in the two ghost layers ----------------------
! X DIRECTION
 count = 1
 DO j = yl, yu
   rho(xug ,j) = rhoG_east_rcv(count)
   rho(xug2,j) = rhoG_east_rcv(count + 1)
   rho(xlg ,j) = rhoG_west_rcv(count)
   rho(xlg2,j) = rhoG_west_rcv(count + 1)
   count = count + 2
 END DO

! Y DIRECTION
 count = 1
 DO i = xl, xu
   rho(i,yug ) = rhoG_north_rcv(count)
   rho(i,yug2) = rhoG_north_rcv(count + 1)
   rho(i,ylg ) = rhoG_south_rcv(count)
   rho(i,ylg2) = rhoG_south_rcv(count + 1)
   count = count + 2
 END DO

! NE DIAGONAL
 rho(xug ,yug ) = rhoG_ne_rcv(1)
 rho(xug2,yug ) = rhoG_ne_rcv(2)
 rho(xug ,yug2) = rhoG_ne_rcv(3)
 rho(xug2,yug2) = rhoG_ne_rcv(4)

! NW DIAGONAL 
 rho(xlg ,yug ) = rhoG_nw_rcv(1)
 rho(xlg2,yug ) = rhoG_nw_rcv(2)
 rho(xlg ,yug2) = rhoG_nw_rcv(3)
 rho(xlg2,yug2) = rhoG_nw_rcv(4)

! SW DIAGONAL 
 rho(xlg ,ylg ) = rhoG_sw_rcv(1)
 rho(xlg2,ylg ) = rhoG_sw_rcv(2)
 rho(xlg ,ylg2) = rhoG_sw_rcv(3)
 rho(xlg2,ylg2) = rhoG_sw_rcv(4)

! SE DIAGONAL
 rho(xug ,ylg ) = rhoG_se_rcv(1)
 rho(xug2,ylg ) = rhoG_se_rcv(2)
 rho(xug ,ylg2) = rhoG_se_rcv(3)
 rho(xug2,ylg2) = rhoG_se_rcv(4)

 RETURN
 END SUBROUTINE MPI_UpdateRhoG

