/*=========================================================================

  Module    : LOCAL STRESS FROM GROMACS TRAJECTORIES
  File      : localstress.h
  Authors   : A. Torres-Sanchez and J. M. Vanegas
  Modified  :
  Purpose   : Compute the local stress from precomputed trajectories in GROMACS
  Date      : 06/23/2015
  Version   :
  Changes   :

     http://www.lacan.upc.edu/LocalStressFromMD

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE. 

     Please, report any bug to either of us:
     juan.m.vanegas@gmail.com
     torres.sanchez.a@gmail.com
=========================================================================*/

// References:
//
// Regarding this program:
// [1] Manual (parent folder/manual)
// [2] J. M. Vanegas; A. Torres-Sanchez; M. Arroyo; J. Chem. Theor. Comput. 10 (2), 691–702 (2014)
// [3] O. H. S. Ollila; H. J. Risselada; M. Louhivuori; E. Lindahl; I. Vattulainen; and S. J. Marrink; Phys. Rev. Lett. 102, 078101 (2009)
//
// General IKN framework and Central Force Decomposition
// [4] E. B. Tadmor; R. E. Miller; Modeling Materials: Continuum, Atomistic and Multiscale Techniques, Cambridge University Press (2011)
// [5] N. C. Admal; E. B. Tadmor; J. Elast. 100, 63 (2010)
//
// Covariant Central Force Decomposition
// [6] A. Torres-Sanchez; J. M. Vanegas; M. Arroyo; Examining the Mechanical Equilibrium of Microscopic Stresses in Molecular Simulations PRL (2015)
// [7] A. Torres-Sanchez; J. M. Vanegas; M. Arroyo; In preparation (2015)
//
// Goetz-Lipowsky Decomposition
// [8] R. Goetz; R. Lipowsky; J. Chem. Phys. 108, 7397 (1998)
//
// Method of planes
// [9] H. Heinz; W. Paul; K. Binder; Phys. Rev. E. 72 066704 (2005)

#ifndef _localstress_h_
#define _localstress_h_


// USEFUL DEFINES
#define enSpat      0
#define enAtom      1


#define encCFD      0
#define enCFD       1
#define enGLD       2
#define enMOP       3
#define enFCD       4
#define enHD_GM     5
#define enHD_LM     6
#define enHD_DIHED  7

#define enSL        0
#define enAll       1
#define enVdw       2
#define enCoul      3
#define enAngles    4
#define enBonds     5
#define enDihp      6
#define enDihi      7
#define enDihrb     8
#define enLincs     9
#define enSettle    10
#define enShake     12
#define enVel       13
#define enNR        14
#define enCMAP      15

#define nRow3       9
#define nCol3       3

#define nRow3_HD    12
#define nCol3_HD    5

#define nRow4       12
#define nCol4       6

#define nRow4_HD    15
#define nCol4_HD    10

#define nRow5       15
#define nCol5       10

#define epsLS       1.0e-10


#include "gmx_lapack.h"
#ifndef F77_FUNC
#define F77_FUNC(name,NAME) name ## _
#endif

// This functions solves an underdetermined or overdetermined system of equations and gives the least squares solution to the problem
void F77_FUNC(dgelsd,DGELSD)( int *, int *, int *, real *, int *, real *, int *, real *, real *, int *, real *, int *, int *, int * );

//MAIN STRUCTURE
typedef struct
{
    int        nAtom;           // Number of atoms
    int        nx;              // Number of grid cells in the x direction
    int        ny;              // Number of grid cells in the y direction
    int        nz;              // Number of grid cells in the z direction
    matrix *   current_grid;    // Grid (either nx*ny*nz or nAtom)
    matrix *   sum_grid;        // Sum Grid
    int        nframes;         // Number of frames
    real       spacing;         // Spacing
    matrix     box;             // Actual box
    matrix     avgbox;          // Average box
    matrix     invbox;          // Inverse of the box
    int        spatatom;        // enSpat or enAtom
    int        contrib;         // which contribution to compute
    int        fdecomp;         // which force decomposition
    gmx_bool   calc_locals;     // true at nstcalcenergy steps */
    gmx_bool   CGlocals;        // Charge-group based local stress */
    int        ePBC;
} gmxLS_locals_grid_t;



//----------------------------------------------------------------------------------------
// gmxLS_distribute_stress
//
// ROOT OF ALL EVIL
// This function reads the gmxLS_locals_grid_t, the number of atoms, the atoms' labels and their
// respective positions and forces, and calls the functions in charge of distributing the stress
// on the grid depending on the local stress flags and the kind of interaction
// Requires:
// grid -> the gmxLS_locals_grid_t where the information concerning the local stress calculation is stored
// nAtom -> number of atoms of the contribution
// atomIDs -> labels of the atoms
// R -> positions of the atoms
// F -> forces on the atoms
void   gmxLS_distribute_stress(gmxLS_locals_grid_t * grid, int nAtom, int* atomIDs, rvec* R, rvec* F);

//----------------------------------------------------------------------------------------
// gmxLS_distribute_interaction
//
// Distributes interactions onto locals_grid (from the initial grid point to the last grid point)
// Requires:
// grid -> the gmxLS_locals_grid_t where the information concerning the local stress calculation is stored
// xi   -> position of particle I (A)
// xj   -> position of particle J (B)
// F    -> pairwise force
// plane -> 0 if all planes contribute, 1 if only the YZ planes contribute,  2 if only the XZ planes contribute,  3 if only the XY planes contribute
void   gmxLS_distribute_interaction(gmxLS_locals_grid_t * grid, rvec xi, rvec xj, rvec F, int plane);

//----------------------------------------------------------------------------------------
// gmxLS_distribute_interaction
//
// Distributes interactions onto locals_grid (from the initial grid point to the last grid point)
// Requires:
// grid    -> the gmxLS_locals_grid_t where the information concerning the local stress calculation is stored
// mass   -> mass of the particle
// atomID -> ID of the atom
// x      -> position of the atom
// va, vb -> velocities from the leap-frog method (we use the average velocity as the one in the middle point)
void gmxLS_spread_vel(gmxLS_locals_grid_t * grid, real mass, int atomID, rvec x, rvec va, rvec vb);

//----------------------------------------------------------------------------------------
// gmxLS_grid_distribute_line_source
//
//Distributes "line sources" onto locals_grid (asuming trilinear weighting functions!)
// This function receives:
// sgrid  -> the grid where we distribute the stress
// a      -> position of particle A
// b      -> position of particle B
// t1     -> starting parametric time (the line segment is parametrized from 0 to 1 from A to B)
// t2     -> ending  parametric time
// x      -> position of the grid point where we are spreading the contribution
// stress -> stress to be spread
// nx, ny, nz -> sizes of the grid in each direction
// gridsp    -> vector with grid spacings in each direction
// invgridsp -> inverse of the volume of the grid cell
// sumfactor -> (this is just for monitoring) part of the line source that has already been spread
// EXPRESSIONS HAVE BEEN COMPUTED ANALYTICALLY!!!
void gmxLS_grid_distribute_line_source(matrix * grid, rvec a, rvec b, real t1, real t2, ivec x, matrix stress, int nx, int ny, int nz, rvec gridsp, real invgridsp, real *sumfactor);

//----------------------------------------------------------------------------------------
// gmxLS_distribute_point_source
//
// Distributes "point sources" onto locals_grid
// Requires:
// grid -> the gmxLS_locals_grid_t where the information concerning the local stress calculation is stored
// pt   -> source point 
// stress -> stress to be distributed
void gmxLS_distribute_point_source(gmxLS_locals_grid_t * grid, rvec pt, matrix stress);

//----------------------------------------------------------------------------------------
//SPECIFIC DECOMPOSITIONS FOR 3,4 AND 5 PARTICLES

// Decompose 3-body potentials (angles)
void   gmxLS_spread_n3(gmxLS_locals_grid_t * grid, rvec Fa, rvec Fb, rvec Fc, rvec Ra, rvec Rb, rvec Rc);

// Decompose Hybrid Decomposition
void   gmxLS_spread_n3_HD(gmxLS_locals_grid_t * grid, rvec Fa, rvec Fb, rvec Fc, rvec Ra, rvec Rb, rvec Rc);

// NOTE: assuming that minimum image convention is employed.
void   gmxLS_get_SDM_min(const rvec r1, const rvec r2, const rvec r3, const rvec r12, const rvec r23, const real sign, rvec ret);
void   gmxLS_get_SDM_gmin(const rvec r1, const rvec r2, const rvec r3, const rvec r12, const rvec r23, rvec ret);
void   gmxLS_get_SDM_lmin(const rvec r1, const rvec r2, const rvec r3, const rvec r12, const rvec r23, rvec ret);
void   gmxLS_get_CM_pos(const rvec r1, const rvec r2, const rvec r3, const rvec r4, rvec ret);
void   gmxLS_get_OFC(const rvec r1, const rvec r2, const rvec dr12, const rvec F2, rvec ret);

// Decompose Settle
void   gmxLS_spread_settle(gmxLS_locals_grid_t * grid, rvec Fa, rvec Fb, rvec Fc, rvec Ra, rvec Rb, rvec Rc);

// Decompose 4-body potentials (dihedrals)
void   gmxLS_spread_n4(gmxLS_locals_grid_t * grid, rvec Fa, rvec Fb, rvec Fc, rvec Fd, rvec Ra, rvec Rb, rvec Rc, rvec Rd);
void gmxLS_spread_n4_HD(gmxLS_locals_grid_t * grid, rvec F1, rvec F2, rvec F3, rvec F4, rvec r1, rvec r2, rvec r3, rvec r4);

// Decompose 5-body potentials (CMAP)
void   gmxLS_spread_n5(gmxLS_locals_grid_t * grid, rvec Fa, rvec Fb, rvec Fc, rvec Fd, rvec Fe, rvec Ra, rvec Rb, rvec Rc, rvec Rd, rvec Re);
//----------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------
// AUXILIARY FUNCTIONS

// Modulo operation
int    gmxLS_modulo (int a, int b);

// Finds the indices on the grid for a given set of coordinates
void   gmxLS_grid_givecoord(gmxLS_locals_grid_t *grid, rvec pt, int *i, int *j, int *k, int pbc);

// Duplicate of grid_givecoord that uses the box directly, needed by g_density3D
void   gmxLS_grid_givecoord2(matrix box, int nx, int ny, int nz, rvec pt, int *i, int *j, int *k, int pbc);

//----------------------------------------------------------------------------------------
// FIVE BODY POTENTIALS -> CALEY-MENGER DERIVATIVES FOR cCFD
// Calculate the derivative of the Caley-Menger determinant for the 5-particles case with respect to d12
double gmxLS_CaleyMenger5Der(double d12,double d13,double d14,double d15,double d23,double d24,double d25,double d34,double d35,double d45);

//Calculate the normal to the shape space for the 5-particles case
void   gmxLS_ShapeSpace5Normal(double d12,double d13,double d14,double d15,double d23,double d24,double d25,double d34,double d35,double d45, double *normal);
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

#endif
