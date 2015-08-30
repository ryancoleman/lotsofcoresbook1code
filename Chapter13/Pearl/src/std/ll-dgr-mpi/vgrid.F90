!-------------------------------------------------------------------------------
! Subroutine : VGrid
! Revision   : 1.1 (2008-11-10)
! Author     : Carlos Rosales Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
!> @file
!! Creates a virtual cartesian grid (2D) linking all MPI-allocated processors.
!> @details
!! Creates a virtual cartesian grid (2D) linking all MPI-allocated processors.
!! The domain limits corresponding to each processor [(xl,xu),(yl,yu)] are
!! are assigned and the rank of each processor and its neigbors in all streaming
!! directions are obtained. By default the processor network is periodic and
!! reordered for faster communication
!! @note
!! Some MPI implementations do not reorder correctly. This may affect the
!! performance but not the results.
!! @warning
!! If you use OpenMPI (v1.2.5) you may have to comment out 'USE MPI' and use
!! instead 'INCLUDE "mpif.h"' because a current bug in OpenMPI prevents the
!! command MPI_CART_CREATE() from being recognized.

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

 SUBROUTINE VGrid

 USE Domain
 USE MPIParams
 USE MPI
 IMPLICIT NONE

! Local variables
 LOGICAL :: reorder
 INTEGER :: complete, direction, partial, shift, mpi_coords(1:mpi_dim)
 INTEGER :: xeast, xwest, ynorth, ysouth
 INTEGER :: MPI_ERR
 INTEGER, DIMENSION(1:mpi_dim) :: dims
 LOGICAL, DIMENSION(1:mpi_dim) :: periodic

! Basic parameters for the partitioning scheme
 dims(1)     = mpi_xdim
 dims(2)     = mpi_ydim
 periodic(1) = .true.
 periodic(2) = .true.
 reorder     = .true.

! Create the new virtual connectivity grid
 CALL MPI_CART_CREATE(MPI_COMM_WORLD, mpi_dim, dims, periodic, reorder, MPI_COMM_VGRID, MPI_ERR)

! Get this processor ID within the virtual grid
 CALL MPI_COMM_RANK(MPI_COMM_VGRID, vproc, MPI_ERR)
 CALL MPI_CART_COORDS(MPI_COMM_VGRID, vproc, mpi_dim, mpi_coords, MPI_ERR)

!--------- Compute domain limits [ (xl,xu) x (yl,yu) ] for this processor ------
! Partitioning in the x direction
 complete = (xmax - xmin) / dims(1)
 partial  = (xmax - xmin) - complete*dims(1)
 IF(mpi_coords(1) + 1 <= partial) THEN
   xl = xmin + (complete + 1)*mpi_coords(1)
   xu = xmin + (complete + 1)*(mpi_coords(1) + 1) - 1
 ELSE
   xl = xmin + complete*mpi_coords(1) + partial
   xu = xmin + complete*(mpi_coords(1) + 1) + partial - 1
 END IF
 IF(MOD(mpi_coords(1) + 1,dims(1)) == 0) xu = xu + 1

! Partitioning in the y direction
 complete = (ymax - ymin) / dims(2)
 partial  = (ymax - ymin) - complete*dims(2)
 IF(mpi_coords(2) + 1 <= partial) THEN
   yl = ymin + (complete + 1)*mpi_coords(2)
   yu = ymin + (complete + 1)*(mpi_coords(2) + 1) - 1
 ELSE
   yl = ymin + complete*mpi_coords(2) + partial
   yu = ymin + complete*(mpi_coords(2) + 1) + partial - 1
 END IF
 IF(MOD(mpi_coords(2) + 1,dims(2)) == 0) yu = yu + 1

! Ghost layer limits (momentum grid)
 xlg  = xl - 1
 xug  = xu + 1 
 ylg  = yl - 1
 yug  = yu + 1

 xlg2 = xl - 2
 xug2 = xu + 2
 ylg2 = yl - 2
 yug2 = yu + 2

! Domain Size (order parameter grid)
 xmin_f = xmin
 ymin_f = ymin
 xmax_f = 2*xmax - 1
 ymax_f = 2*ymax - 1

 xl_f = 2*xl - 2
 xu_f = 2*xu - 1
 yl_f = 2*yl - 2
 yu_f = 2*yu - 1

! Ghost layer limits (order parameter grid)
 xlg_f  = xl_f - 1
 xug_f  = xu_f + 1
 ylg_f  = yl_f - 1
 yug_f  = yu_f + 1

 xlg2_f = xl_f - 2
 xug2_f = xu_f + 2
 ylg2_f = yl_f - 2
 yug2_f = yu_f + 2

 xlg3_f = xl_f - 3
 xug3_f = xu_f + 3
 ylg3_f = yl_f - 3
 yug3_f = yu_f + 3

!--------- Determine neighbours of this processor ------------------------------
! MPI_CART counts dimensions using 0-based arithmetic so that
! direction = 0 -> x  |  direction = 1 -> y
 shift     = 1
 direction = 0
 CALL MPI_CART_SHIFT(MPI_COMM_VGRID, direction, shift, west, east, MPI_ERR)
 direction = 1
 CALL MPI_CART_SHIFT(MPI_COMM_VGRID, direction, shift, south, north, MPI_ERR)

! Obtain coordinates of neighbours of this processor in the x and y directions
 CALL MPI_CART_COORDS(MPI_COMM_VGRID, east, mpi_dim, mpi_coords, MPI_ERR)
 xeast = mpi_coords(1)
 CALL MPI_CART_COORDS(MPI_COMM_VGRID, west, mpi_dim, mpi_coords, MPI_ERR)
 xwest = mpi_coords(1)
 CALL MPI_CART_COORDS(MPI_COMM_VGRID, north, mpi_dim, mpi_coords, MPI_ERR)
 ynorth = mpi_coords(2)
 CALL MPI_CART_COORDS(MPI_COMM_VGRID, south, mpi_dim, mpi_coords, MPI_ERR)
 ysouth = mpi_coords(2)

! Obtain the ranks of the diagonal neighbours
 mpi_coords(1) = xeast
 mpi_coords(2) = ynorth
 CALL MPI_CART_RANK(MPI_COMM_VGRID,mpi_coords,ne,MPI_ERR)
 mpi_coords(1) = xwest
 mpi_coords(2) = ysouth
 CALL MPI_CART_RANK(MPI_COMM_VGRID,mpi_coords,sw,MPI_ERR)
 mpi_coords(1) = xwest
 mpi_coords(2) = ynorth
 CALL MPI_CART_RANK(MPI_COMM_VGRID,mpi_coords,nw,MPI_ERR)
 mpi_coords(1) = xeast
 mpi_coords(2) = ysouth
 CALL MPI_CART_RANK(MPI_COMM_VGRID,mpi_coords,se,MPI_ERR)

 RETURN
 END SUBROUTINE VGrid
