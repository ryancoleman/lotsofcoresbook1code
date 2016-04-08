!-------------------------------------------------------------------------------
! Subroutine : VGrid
! Revision   : 1.3 (2013/11/12)
! Author     : Carlos Rosales-Fernandez [carlos.rosales.fernandez(at)gmail.com]
!-------------------------------------------------------------------------------
! Creates a virtual cartesian grid (3D) linking all MPI-allocated processors.
! The domain limits corresponding to each processor [(xl,xu),(yl,yu),(zl, zu)]
! are assigned and the rank of each processor and its neigbors in all streaming
! directions are obtained. By default the processor network is periodic and
! reordered for faster communication
!
! *** NOTE ***
! Many MPI implementations do not reorder correctly. This may affect the
! performance but not the results.
!
! *** WARNING ***
! If you use old OpenMPI versions you may have to comment out 'USE MPI' and use
! instead 'INCLUDE "mpif.h"' because an old bug in OpenMPI prevents the
! command MPI_CART_CREATE() from being recognized.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Copyright 2013 Carlos Rosales Fernandez and The University of Texas at Austin.
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

!  Common Variables
 USE Domain
 USE MPIParams
 USE MPI
 IMPLICIT NONE

! Local variables
 INTEGER :: complete, direction, partial, shift
 INTEGER :: MPI_ERR, x, xeast, xwest, y, ynorth, ysouth, z, ztop, zbottom
 INTEGER, DIMENSION(1:mpi_dim) :: dims, mpi_coords
 LOGICAL, DIMENSION(1:mpi_dim) :: periodic
 LOGICAL :: reorder

! Initialize data for domain partitioning. Defaults are:
! Partitioning is periodic in all dimensions (periodic = .true.)
! CPUs are reordered in the grid for proximity (reorder = .true.)
 dims(1)     = mpi_xdim
 dims(2)     = mpi_ydim
 dims(3)     = mpi_zdim
 periodic(1) = .true.
 periodic(2) = .true.
 periodic(3) = .true.
 reorder     = .true.

! Create the new virtual connectivity grid
 CALL MPI_CART_CREATE(MPI_COMM_WORLD, mpi_dim, dims, periodic, reorder, MPI_COMM_VGRID, MPI_ERR)

! Get this processor ID and coordinates within the virtual grid
 CALL MPI_COMM_RANK(MPI_COMM_VGRID, vproc, MPI_ERR)
 CALL MPI_CART_COORDS(MPI_COMM_VGRID, vproc, mpi_dim, mpi_coords, MPI_ERR)
 x = mpi_coords(1)
 y = mpi_coords(2)
 z = mpi_coords(3)

!-------- Compute limits [ (xl,xu) x (yl,yu) x (zl,zu)] for this processor -----
!  Partitioning in the x direction
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
!  Partitioning in the y direction
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
!  Partitioning in the y direction
 complete = (zmax - zmin) / dims(3)
 partial  = (zmax - zmin) - complete*dims(3)
 IF(mpi_coords(3) + 1 <= partial) THEN
   zl = zmin + (complete + 1)*mpi_coords(3)
   zu = zmin + (complete + 1)*(mpi_coords(3) + 1) - 1
 ELSE
   zl = zmin + complete*mpi_coords(3) + partial
   zu = zmin + complete*(mpi_coords(3) + 1) + partial - 1
 END IF
 IF(MOD(mpi_coords(3) + 1,dims(3)) == 0) zu = zu + 1

! Modify limits to include ghost layers
 xlg = xl - 1
 xug = xu + 1
 ylg = yl - 1
 yug = yu + 1
 zlg = zl - 1
 zug = zu + 1

!-------- Determine neighbours of this processor -------------------------------
! MPI_CART counts dimensions using 0-based arithmetic so that
! direction = 0 -> x  |  direction = 1 -> y  |  direction = 2 -> z
 shift     = 1
 direction = 0
 CALL MPI_CART_SHIFT(MPI_COMM_VGRID, direction, shift, west, east, MPI_ERR)
 direction = 1
 CALL MPI_CART_SHIFT(MPI_COMM_VGRID, direction, shift, south, north, MPI_ERR)
 direction = 2
 CALL MPI_CART_SHIFT(MPI_COMM_VGRID, direction, shift, bot, top, MPI_ERR)

!  Coordinates of the neighbours of this processor in the x, y, and z directions
 CALL MPI_CART_COORDS(MPI_COMM_VGRID, east, mpi_dim, mpi_coords, MPI_ERR)
 xeast = mpi_coords(1)
 CALL MPI_CART_COORDS(MPI_COMM_VGRID, west, mpi_dim, mpi_coords, MPI_ERR)
 xwest = mpi_coords(1)
 CALL MPI_CART_COORDS(MPI_COMM_VGRID, north, mpi_dim, mpi_coords, MPI_ERR)
 ynorth = mpi_coords(2)
 CALL MPI_CART_COORDS(MPI_COMM_VGRID, south, mpi_dim, mpi_coords, MPI_ERR)
 ysouth = mpi_coords(2)
 CALL MPI_CART_COORDS(MPI_COMM_VGRID, top, mpi_dim, mpi_coords, MPI_ERR)
 ztop = mpi_coords(3)
 CALL MPI_CART_COORDS(MPI_COMM_VGRID, bot, mpi_dim, mpi_coords, MPI_ERR)
 zbottom = mpi_coords(3)

! Get the ranks of the diagonal neighbours
!  Neighbors along x edges
 mpi_coords(1) = x
 mpi_coords(2) = ynorth
 mpi_coords(3) = ztop
 CALL MPI_CART_RANK(MPI_COMM_VGRID,mpi_coords,tn,MPI_ERR)
 mpi_coords(1) = x
 mpi_coords(2) = ysouth
 mpi_coords(3) = ztop
 CALL MPI_CART_RANK(MPI_COMM_VGRID,mpi_coords,ts,MPI_ERR)
 mpi_coords(1) = x
 mpi_coords(2) = ysouth
 mpi_coords(3) = zbottom
 CALL MPI_CART_RANK(MPI_COMM_VGRID,mpi_coords,bs,MPI_ERR)
 mpi_coords(1) = x
 mpi_coords(2) = ynorth
 mpi_coords(3) = zbottom
 CALL MPI_CART_RANK(MPI_COMM_VGRID,mpi_coords,bn,MPI_ERR)

!  Neighbors along y edges
 mpi_coords(1) = xeast
 mpi_coords(2) = y
 mpi_coords(3) = ztop
 CALL MPI_CART_RANK(MPI_COMM_VGRID,mpi_coords,te,MPI_ERR)
 mpi_coords(1) = xwest
 mpi_coords(2) = y
 mpi_coords(3) = ztop
 CALL MPI_CART_RANK(MPI_COMM_VGRID,mpi_coords,tw,MPI_ERR)
 mpi_coords(1) = xwest
 mpi_coords(2) = y
 mpi_coords(3) = zbottom
 CALL MPI_CART_RANK(MPI_COMM_VGRID,mpi_coords,bw,MPI_ERR)
 mpi_coords(1) = xeast
 mpi_coords(2) = y
 mpi_coords(3) = zbottom
 CALL MPI_CART_RANK(MPI_COMM_VGRID,mpi_coords,be,MPI_ERR)

!  Neighbors along z edges
 mpi_coords(1) = xeast
 mpi_coords(2) = ynorth
 mpi_coords(3) = z
 CALL MPI_CART_RANK(MPI_COMM_VGRID,mpi_coords,ne,MPI_ERR)
 mpi_coords(1) = xwest
 mpi_coords(2) = ysouth
 mpi_coords(3) = z
 CALL MPI_CART_RANK(MPI_COMM_VGRID,mpi_coords,sw,MPI_ERR)
 mpi_coords(1) = xwest
 mpi_coords(2) = ynorth
 mpi_coords(3) = z
 CALL MPI_CART_RANK(MPI_COMM_VGRID,mpi_coords,nw,MPI_ERR)
 mpi_coords(1) = xeast
 mpi_coords(2) = ysouth
 mpi_coords(3) = z
 CALL MPI_CART_RANK(MPI_COMM_VGRID,mpi_coords,se,MPI_ERR)

! Set domain limits
 NX = xu - xl + 1
 NY = yu - yl + 1
 NZ = zu - zl + 1

! Set ghost limits
 NXG = NX + 2
 NYG = NY + 2
 NZG = NZ + 2
 NG  = NXG*NYG*NZG

 RETURN
 END SUBROUTINE VGrid
