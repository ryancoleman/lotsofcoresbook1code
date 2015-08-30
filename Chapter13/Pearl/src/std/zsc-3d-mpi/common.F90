!-------------------------------------------------------------------------------
! File     : Common
! Revision : 1.0 (2008-06-15)
! Author   : David S. Whyte [david(at)ihpc.a-star.edu.sg]
!-------------------------------------------------------------------------------
!> @file
!! Modules that contain common variables

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

!> @brief Definition of single and double data types
 MODULE NTypes
 IMPLICIT NONE
 SAVE

 INTEGER, PARAMETER :: SGL = KIND(1.0)
 INTEGER, PARAMETER :: DBL = KIND(1.D0)

 END MODULE NTypes

!> @brief Parameters related to the geometry and the time intervals
 MODULE Domain
 USE NTypes, ONLY : DBL
 IMPLICIT NONE
 SAVE

! Time step and limits
 INTEGER :: iStep, tCall, tDump, tStat, MaxStep, RelaxStep

! Domain size
! (xmax, xmin) : computational domain size limits
! (xl  , xu  ) : local minimun and maximum values for x in a given processor
! (xlg , xug ) : local ghost node position (xlg = xl-1, xug = xu+1)
 INTEGER :: xjump, xo, yo, zo
 INTEGER :: now, nxt, NX, NY, NZ, NPX, NPy, NPZ
 INTEGER :: xmin, xmax, xl, xu, xlg, xug
 INTEGER :: ymin, ymax, yl, yu, ylg, yug
 INTEGER :: zmin, zmax, zl, zu, zlg, zug
 INTEGER :: xsize, xsize2, xsize3, xsize4, xsize5
 INTEGER :: ysize, ysize2, ysize3, ysize4, ysize5
 INTEGER :: zsize, zsize2, zsize3, zsize4, zsize5
 INTEGER :: xedge, yedge, zedge

! ndim = spatial dimension
 INTEGER, PARAMETER :: ndim = 3

! Array for neigbors
 INTEGER, ALLOCATABLE, DIMENSION(:,:,:,:) :: ni

! Inverse of the initial volume of the discrete phase 
 REAL(KIND = DBL) :: invInitVol

! Define constants used in the code (inv6 = 1/6, inv12 = 1/12, invPi = 1/pi)
 REAL(KIND = DBL), PARAMETER :: inv6  = 0.16666666666666666667D0
 REAL(KIND = DBL), PARAMETER :: inv12 = 0.08333333333333333333D0
 REAL(KIND = DBL), PARAMETER :: invPi = 0.31830988618379067154D0

 END MODULE Domain

!> @brief Parameters related to hydrodynamics quantities
 MODULE FluidParams
 USE NTypes, ONLY : DBL
 IMPLICIT NONE
 SAVE

 INTEGER :: nBubbles
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:,:) :: bubbles
 REAL(KIND = DBL), DIMENSION(1:11) :: Convergence
 REAL(KIND = DBL) :: sigma, IntWidth, Gamma
 REAL(KIND = DBL) :: alpha, alpha4, kappa
 REAL(KIND = DBL) :: tauPhi, invTauPhi, invTauPhi1, phiStar, phiStar2, phiStar4
 REAL(KIND = DBL) :: rhoL, rhoH
 REAL(KIND = DBL) :: tauRho, invTauRho, invTauRhoOne, invTauRhoHalf
 REAL(KIND = DBL) :: eta, eta2, invEta2
 REAL(KIND = DBL) :: gravity, grav, Eo, eps, pConv

 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:,:,:) :: phi

 END MODULE FluidParams

!> @brief Parameters related to the LBM discretization
 MODULE LBMParams
 USE NTypes, ONLY : DBL
 IMPLICIT NONE
 SAVE

! fdim = order parameter dimension - 1 (D3Q7)
! gdim = momentum dimension - 1        (D3Q19)
 INTEGER, PARAMETER :: fdim = 6
 INTEGER, PARAMETER :: gdim = 18

! Distribution functions
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: f, g

! D3Q7/D3Q19 Lattice speed of sound (Cs = 1/DSQRT(3), Cs_sq = 1/3, invCs_sq = 3)
 REAL(KIND = DBL), PARAMETER :: Cs       = 0.57735026918962576451D0
 REAL(KIND = DBL), PARAMETER :: Cs_sq    = 0.33333333333333333333D0
 REAL(KIND = DBL), PARAMETER :: invCs_sq = 3.00000000000000000000D0

! Distributions weights (D3Q19: Eg0 = 1/3, Eg1 = 1/18, Eg2 = 1/36)
 REAL(KIND = DBL), PARAMETER :: Eg0 = 0.33333333333333333333D0
 REAL(KIND = DBL), PARAMETER :: Eg1 = 0.05555555555555555556D0
 REAL(KIND = DBL), PARAMETER :: Eg2 = 0.02777777777777777778D0

! Modified distribution weight (to avoid operations: EgiC = Egi*invCs_sq)
 REAL(KIND = DBL), PARAMETER :: Eg0C = 1.00000000000000000000D0
 REAL(KIND = DBL), PARAMETER :: Eg1C = 0.16666666666666666667D0
 REAL(KIND = DBL), PARAMETER :: Eg2C = 0.08333333333333333333D0
 REAL(KIND = DBL) :: w, Eg0n, Eg1n, Eg2n, EgC0n, EgC1n, EgC2n

 END MODULE LBMParams

!> @brief MPI related parameters and information exchange buffer arrays
!> @details
!! @param nprocs  : total number of processors
!! @param proc    : processor ID
!! @param vproc   : processor ID in virtual grid
!! @param mpi_dim : virtual grid partition scheme (1->stripes, 2->boxes, 3->cubes)
 MODULE MPIParams
 USE NTypes, ONLY : DBL
 IMPLICIT NONE
 SAVE

! Communication parameters
 INTEGER :: nprocs, proc, vproc, MPI_COMM_VGRID
 INTEGER :: mpi_xdim, mpi_ydim, mpi_zdim
 INTEGER :: east, west, north, south, top, bot
 INTEGER :: ne, nw, se, sw, te, tw, be, bw, tn, ts, bn, bs
 INTEGER, PARAMETER :: master  = 0
 INTEGER, PARAMETER :: mpi_dim = 3
 INTEGER, PARAMETER :: TAG1 = 1, TAG2 = 2, TAG3 = 3, TAG4 = 4

! Information exchange buffers (x direction)
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: f_west_snd, f_east_snd
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: f_west_rcv, f_east_rcv
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: g_west_snd, g_east_snd
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: g_west_rcv, g_east_rcv
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: phi_west_snd, phi_east_snd
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: phi_west_rcv, phi_east_rcv

! Information exchange buffers (y direction)
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: f_south_snd, f_north_snd
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: f_south_rcv, f_north_rcv
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: g_south_snd, g_north_snd
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: g_south_rcv, g_north_rcv
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: phi_south_snd, phi_north_snd
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: phi_south_rcv, phi_north_rcv

! Information exchange buffers (z direction)
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: f_bot_snd, f_top_snd
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: f_bot_rcv, f_top_rcv
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: g_bot_snd, g_top_snd
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: g_bot_rcv, g_top_rcv
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: phi_bot_snd, phi_top_snd
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: phi_bot_rcv, phi_top_rcv

! Information exchange buffers (diagonals x)
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: g_tn_snd, g_tn_rcv
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: g_bn_snd, g_bn_rcv
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: g_ts_snd, g_ts_rcv
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: g_bs_snd, g_bs_rcv
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: phi_tn_snd, phi_tn_rcv
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: phi_bn_snd, phi_bn_rcv
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: phi_ts_snd, phi_ts_rcv
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: phi_bs_snd, phi_bs_rcv

! Information exchange buffers (diagonals y)
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: g_te_snd, g_te_rcv
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: g_be_snd, g_be_rcv
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: g_tw_snd, g_tw_rcv
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: g_bw_snd, g_bw_rcv
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: phi_te_snd, phi_te_rcv
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: phi_be_snd, phi_be_rcv
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: phi_tw_snd, phi_tw_rcv
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: phi_bw_snd, phi_bw_rcv

! Information exchange buffers (diagonals z)
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: g_ne_snd, g_ne_rcv
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: g_se_snd, g_se_rcv
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: g_nw_snd, g_nw_rcv
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: g_sw_snd, g_sw_rcv
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: phi_ne_snd, phi_ne_rcv
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: phi_se_snd, phi_se_rcv
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: phi_nw_snd, phi_nw_rcv
 REAL(KIND = DBL), ALLOCATABLE, DIMENSION(:) :: phi_sw_snd, phi_sw_rcv

 END MODULE MPIParams
