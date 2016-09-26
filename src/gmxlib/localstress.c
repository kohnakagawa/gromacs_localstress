/*=========================================================================

  Module    : LOCAL STRESS FROM GROMACS TRAJECTORIES
  File      : localstress.c
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
// [6] A. Torres-Sanchez; J. M. Vanegas; M. Arroyo; Examining the Mechanical Equilibrium of Microscopic Stresses in Molecular Simulations (2015)
// [7] A. Torres-Sanchez; J. M. Vanegas; M. Arroyo; In preparation (2015)
//
// Goetz-Lipowsky Decomposition
// [8] R. Goetz; R. Lipowsky; J. Chem. Phys. 108, 7397 (1998)
//
// Method of planes
// [9] H. Heinz; W. Paul; K. Binder; Phys. Rev. E. 72 066704 (2005)

#include "localstress.h"
#include "types/simple.h"
#include "typedefs.h"
#include "vec.h"
#include "pbc.h"
#include <stdio.h>
#include <math.h>
#include <assert.h>

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
void gmxLS_distribute_stress(gmxLS_locals_grid_t * grid, int nAtom, int* atomIDs, rvec* R, rvec* F)
{
    int nDim = 3;
    int n;
    int i,j;
    real temp;

    matrix *sgrid;
    matrix stress;

    if(grid==NULL || grid->calc_locals==FALSE)
        return;

    // If spatatom==enSpat distribute the stress spatially following Noll's procedure
    if (grid->spatatom == enSpat)
    {
        // Depending on the number of atoms, call a different gmxLS_spread function
        switch (nAtom)
        {

            case 2:
                gmxLS_distribute_interaction(grid, R[0], R[1], F[0],0);
                break;

            case 3:
                gmxLS_spread_n3(grid, F[0], F[1], F[2], R[0], R[1], R[2]);
                break;

            case -3:
                gmxLS_spread_settle(grid, F[0], F[1], F[2], R[0], R[1], R[2]);
                break;

            case 4:
                gmxLS_spread_n4(grid, F[0], F[1], F[2], F[3], R[0], R[1], R[2], R[3]);
                break;

            case 5:
                gmxLS_spread_n5(grid, F[0], F[1], F[2], F[3], F[4], R[0], R[1], R[2], R[3], R[4]);
                break;

            // NOTE: To compare with the original implementation implemented
            //       by Torres-Sanchez, et al., we prepare special case.
            //       6 does not correspond to the number of atoms.
            case 6:
                gmxLS_spread_n3_HD(grid, F[0], F[1], F[2], R[0], R[1], R[2]);
                break;

            default:
                break;
        }

    }
    // If spatatom==enAtom, distributes the stress per atom using the atomic stress definition
    else if (grid->spatatom == enAtom)
    {
        // This is because SETTLE calls the function with nAtom=-3
        if (nAtom < 0) nAtom = -nAtom;

        for (n = 0; n < nAtom; n++ )
        {
            if ( atomIDs[n] >= grid->nAtom )
            {
                printf("ERROR:: the atom label %d is equal or larger than the total number of atoms %d\n",atomIDs[n], grid->nAtom);
                return;
            }
        }

        sgrid = grid->current_grid;

        //Initialize the value of the (local) stress to 0
        for( i = 0; i< nDim; i++ )
        {
            stress[i][i] = 0.0;
            for( j=i+1; j< nDim; j++ )
            {
                stress[i][j] = 0.0;
                stress[j][i] = 0.0;
            }
        }


        //Calculate the stress
        for( n = 0; n < nAtom; n++ )
        {
            for( i = 0; i< nDim; i++ )
            {
                temp = -(F[n][i] * R[n][i])/nAtom;
                stress[i][i] += temp;
                for( j=i+1; j< nDim; j++ )
                {
                    temp = -(F[n][i] * R[n][j])/nAtom;
                    stress[i][j] += temp;
                    stress[j][i] += temp;
                }
            }
        }

        for ( n = 0; n < nAtom; n++ )
        {
            m_add(sgrid[atomIDs[n]], stress, sgrid[atomIDs[n]]);
        }
    }

    return;
}

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
// the plane variable is used for the stress on geometric centers option, see ref. [9]
void gmxLS_distribute_interaction(gmxLS_locals_grid_t * grid, rvec xi, rvec xj, rvec F, int plane)
{
    int nDim = 3;

    real gridp, oldt, sumfactor, invgridsp;

    rvec box, gridsp, curpt, t, d_cgrid, diff;

    int  i, j, k, ii, jj, kk; //counters
    int nx, ny, nz;           //size of the box
    int tof;                  //true or false

    ivec i1; //grid cell corresponding to particle I (A)
    ivec i2; //grid cell corresponding to particle J (B)
    ivec x;  //cell during spreding
    ivec xn; //next cell during spreading
    ivec c; //director

    matrix gridres, *sgrid, stress;

    int timesInLoop; //this is to avoid that the program gets stacked in an infinite loop

    // Load grid parameters and calculate spacing
    sgrid = grid->current_grid;
    box[XX] = grid->box[XX][XX];
    box[YY] = grid->box[YY][YY];
    box[ZZ] = grid->box[ZZ][ZZ];
    nx = grid->nx;
    ny = grid->ny;
    nz = grid->nz;
    gridsp[XX] = box[XX]/nx;
    gridsp[YY] = box[YY]/ny;
    gridsp[ZZ] = box[ZZ]/nz;
    invgridsp = 1/(gridsp[0]*gridsp[1]*gridsp[2]);

    //------------------------------------------------------------------------------------
    // Calculate the stress tensor

    // difference between xj and xi
    rvec_sub(xj,xi,diff);

    // The use of planes is only necessary for the stress on geometric centers option see ref. [9]
    if( plane > 0 ) //If not all planes contribute, set all to 0 to avoid corrupted data
    {
        stress[XX][XX] = 0.0;
        stress[YY][XX] = 0.0;
        stress[ZZ][XX] = 0.0;
        stress[XX][YY] = 0.0;
        stress[YY][YY] = 0.0;
        stress[ZZ][YY] = 0.0;
        stress[XX][ZZ] = 0.0;
        stress[YY][ZZ] = 0.0;
        stress[ZZ][ZZ] = 0.0;
    }
    if( plane == 1 || plane==0 )
    {
        stress[XX][XX] = F[0]*diff[XX];
        stress[YY][XX] = F[1]*diff[XX];
        stress[ZZ][XX] = F[2]*diff[XX];
    }
    if( plane == 2 || plane==0 )
    {
        stress[XX][YY] = F[0]*diff[YY];
        stress[YY][YY] = F[1]*diff[YY];
        stress[ZZ][YY] = F[2]*diff[YY];
    }
    if( plane == 3 || plane==0 )
    {
        stress[XX][ZZ] = F[0]*diff[ZZ];
        stress[YY][ZZ] = F[1]*diff[ZZ];
        stress[ZZ][ZZ] = F[2]*diff[ZZ];
    }
    //------------------------------------------------------------------------------------

    //------------------------------------------------------------------------------------
    //Distribute the stress

    // Calculate the grid coordinates (no pbc) for the extreme points
    gmxLS_grid_givecoord(grid,xi, &i1[0], &i1[1], &i1[2], 0);
    gmxLS_grid_givecoord(grid,xj, &i2[0], &i2[1], &i2[2], 0);

    // d_cgrid = vector from the center of the present cell to the initial point
    // c is a vector that guide the advance in each coordinate (+1 if it has to advance in this coordinate, -1 if it has to go back or 0 if it has to do nothing)
    for(i = 0; i < nDim; i++) 
    {
        d_cgrid[i] = xi[i]-(i1[i]+0.5)*gridsp[i];
        if(i2[i]>i1[i])
        {
            c[i] = 1;
        } 
        else if(i1[i]>i2[i])
        {
            c[i]=-1;
        } else
        {
            c[i] = 0;
        }
    }

    // First cross point with aplane (if there is at least one i.e. c[i] != 0)
    for(i = 0; i < nDim; i++ )
    {
        if(c[i]==1)
        {
            xn[i] = i1[i]+1;                    //label of the next cell is 1 step further than the previous in this direction
            gridp = xn[i] * gridsp[i];          //position of the cell
            t[i] = (xi[i]-gridp)/(xi[i]-xj[i]); //parametric time of crossing
        } 
        else if(c[i]==-1)
        {
            xn[i] = i1[i];                      //label of the next cell is the same in this direction
            gridp = xn[i] * gridsp[i];          //position of the cell
            t[i] = (xi[i]-gridp)/(xi[i]-xj[i]); //parametric time of crossing
        } 
        else 
            t[i]=1.1;                           //This sets the time larger than 1 because there is no crossing

        x[i] = i1[i];

    }

    oldt = 0;      //previous parametric time of crossing
    sumfactor = 0; //This is just to check that the contribution has been completely spread

    // While we don't reach the last point...

    tof = (c[0]*x[0]<c[0]*i2[0]+epsLS)||(c[1]*x[1]<c[1]*i2[1]+epsLS)||(c[2]*x[2]<c[2]*i2[2]+epsLS);

    timesInLoop = 0;

    while(tof)
    {

        if(t[0]<t[1]+epsLS && t[0]<t[2]+epsLS)
        {
            // Distribute the contribution
            gmxLS_grid_distribute_line_source(sgrid,diff,d_cgrid,oldt,t[0],x,stress,nx,ny,nz,gridsp,invgridsp,&sumfactor);

            // Move
            d_cgrid[0] -= c[0] * gridsp[0];
            oldt = t[0];
            
            x[0] += c[0];
            xn[0] += c[0];

            // Next cross point:
            gridp = xn[0] *gridsp[0];
            t[0] = (xi[0]-gridp)/(xi[0]-xj[0]);
        } 
        else if (t[1]<t[0]+epsLS && t[1]<t[2]+epsLS)
        {
            // Distribute the contribution
            gmxLS_grid_distribute_line_source(sgrid,diff,d_cgrid,oldt,t[1],x,stress,nx,ny,nz,gridsp,invgridsp,&sumfactor);

            // Move
            d_cgrid[1] -= c[1] * gridsp[1];
            oldt = t[1];

            x[1]+=c[1];
            xn[1]+=c[1];

            // Next cross point:
            gridp = xn[1] * gridsp[1];
            t[1] = (xi[1]-gridp)/(xi[1]-xj[1]);

        } 
        else if (t[2]<t[0]+epsLS && t[2]<t[1]+epsLS)
        {
            // Distribute the contribution
            gmxLS_grid_distribute_line_source(sgrid,diff,d_cgrid,oldt,t[2],x,stress,nx,ny,nz,gridsp,invgridsp,&sumfactor);

            // Move
            d_cgrid[2] -= c[2] * gridsp[2];
            oldt  = t[2];

            x[2]+=c[2];
            xn[2]+=c[2];

            // Next cross point:
            gridp = xn[2] * gridsp[2];
            t[2] = (xi[2]-gridp)/(xi[2]-xj[2]);
        }
        else
        {
            printf("ERROR::gmx_spread_local_stress_on_grid: t=(%lf, %lf, %lf), i1=(%d, %d, %d), i2=(%d, %d, %d), x=(%d, %d, %d)\n",t[0], t[1], t[2], i1[0], i1[1], i1[2], i2[0], i2[1], i2[2], x[0], x[1], x[2]);
            return;
        }

        tof = (c[0]*x[0]<c[0]*i2[0])||(c[1]*x[1]<c[1]*i2[1])||(c[2]*x[2]<c[2]*i2[2]);
        timesInLoop ++;
        
        if(timesInLoop > 10000000) //To avoid infinite loops
        {
            return;
        }
    }

    // Distribute the last contribution

    gmxLS_grid_distribute_line_source(sgrid,diff,d_cgrid,oldt,1,x,stress,nx,ny,nz,gridsp,invgridsp,&sumfactor);
    //------------------------------------------------------------------------------------
}

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
void gmxLS_spread_vel(gmxLS_locals_grid_t * grid, real mass, int atomID, rvec x, rvec va, rvec vb)
{

    matrix *sgrid;
    matrix stress;

    stress[XX][XX] = -mass*(va[XX]*va[XX] + vb[XX]*vb[XX])/2;
    stress[XX][YY] = -mass*(va[XX]*va[YY] + vb[XX]*vb[YY])/2;
    stress[XX][ZZ] = -mass*(va[XX]*va[ZZ] + vb[XX]*vb[ZZ])/2;
    stress[YY][XX] = -mass*(va[YY]*va[XX] + vb[YY]*vb[XX])/2;
    stress[YY][YY] = -mass*(va[YY]*va[YY] + vb[YY]*vb[YY])/2;
    stress[YY][ZZ] = -mass*(va[YY]*va[ZZ] + vb[YY]*vb[ZZ])/2;
    stress[ZZ][XX] = -mass*(va[ZZ]*va[XX] + vb[ZZ]*vb[XX])/2;
    stress[ZZ][YY] = -mass*(va[ZZ]*va[YY] + vb[ZZ]*vb[YY])/2;
    stress[ZZ][ZZ] = -mass*(va[ZZ]*va[ZZ] + vb[ZZ]*vb[ZZ])/2;

    // If spatatom==enSpat distribute the stress spatially following Noll's procedure
    if (grid->spatatom == enSpat)
    {
        gmxLS_distribute_point_source(grid, x, stress);
    }
    else if (grid->spatatom == enAtom)
    {
        sgrid = grid->current_grid;
        m_add(sgrid[atomID],stress,sgrid[atomID]);
    }
}

//----------------------------------------------------------------------------------------
// gmxLS_grid_distribute_line_source
//
// Distributes "line sources" onto a grid point (asuming trilinear weighting functions!) 
// Requires:
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
void gmxLS_grid_distribute_line_source(matrix * sgrid, rvec a, rvec b, real t1, real t2, ivec x, matrix stress, int nx, int ny, int nz, rvec gridsp, real invgridsp, real *sumfactor)
{

    real factor;
    real dummy1,dummy2,dummy3,dummy4,dummy5,dummy6,dummy7,dummy8,dummy9,dummy10,dummy11,dummy12;
    int i,j,k;
    int ii,jj,kk;
    matrix gridres;
    real t12,t22,t13,t23;

    t12=t1*t1;
    t13=t12*t1;
    t22=t2*t2;
    t23=t22*t2;

    ii=x[0]; jj=x[1]; kk=x[2];

    dummy1 = -2*a[0]*a[1]*a[2]*t12*t12;
    dummy2 =  2*a[0]*a[1]*a[2]*t22*t22;

    for(i=1;i>=-1;i-=2)
    {
        ii+=i;
        dummy3 = 2*b[0]+i*gridsp[0];
        dummy7 = a[1]*a[2]*dummy3;

        for(j=1;j>=-1;j-=2)
        {
            jj+=j;
            dummy4 = 2*b[1]+j*gridsp[1];
            dummy6 = dummy3*dummy4;
            dummy8 = a[0]*a[2]*dummy4;
            dummy10= a[2]*dummy3*dummy4;

            for(k=1;k>=-1;k-=2)
            {
                kk+=k;
                dummy5 = 2*b[2]+k*gridsp[2];
                dummy9 = a[1]*a[0]*dummy5;
                dummy11= a[1]*dummy3*dummy5;
                dummy12= a[0]*dummy4*dummy5;
                factor = i*j*k*0.125*invgridsp*invgridsp*(dummy1+dummy2+(t2-t1)*dummy6*dummy5+1.333333333333*(t23-t13)
                            *(dummy7+dummy8+dummy9)+(t22-t12)*(dummy10+dummy11+dummy12));

                *sumfactor=*sumfactor+factor;
                msmul(stress,factor,gridres);
                m_add(sgrid[gmxLS_modulo(ii,nx)*nz*ny+gmxLS_modulo(jj,ny)*nz+gmxLS_modulo(kk,nz)], gridres,
                        sgrid[gmxLS_modulo(ii,nx)*nz*ny+gmxLS_modulo(jj,ny)*nz+gmxLS_modulo(kk,nz)]);
            }
        }
    }

}


//----------------------------------------------------------------------------------------
// gmxLS_distribute_point_source
//
// Distributes "point sources" onto locals_grid
// Requires:
// grid -> the gmxLS_locals_grid_t where the information concerning the local stress calculation is stored
// pt   -> source point 
// stress -> stress to be distributed
void gmxLS_distribute_point_source(gmxLS_locals_grid_t * grid, rvec pt, matrix stress)
{

    // Spreads the velocity in one point

    rvec box, gridsp;
    int i, j, k, ii, jj, kk, iii, jjj, kkk, nx, ny, nz;
    matrix gridres, *sgrid;
    real factor, invgridsp, dummy1, dummy2;

    // Load parameters
    sgrid = grid->current_grid;
    nx = grid->nx;
    ny = grid->ny;
    nz = grid->nz;
    box[XX] = grid->box[XX][XX];
    box[YY] = grid->box[YY][YY];
    box[ZZ] = grid->box[ZZ][ZZ];

    // Define the grid spacing in each direction
    gridsp[XX] = box[XX]/nx;
    gridsp[YY] = box[YY]/ny;
    gridsp[ZZ] = box[ZZ]/nz;
    invgridsp = 1/(gridsp[XX]*gridsp[YY]*gridsp[ZZ]);

    // Get the coordinates of the point in the grid
    gmxLS_grid_givecoord(grid,pt,&ii,&jj,&kk,0);
    iii=ii;
    jjj=jj;
    kkk=kk;

    // Spread it
    for(i=1;i>=-1;i-=2)
    {
        iii+=i;
        dummy1 = i * invgridsp * invgridsp * (pt[0]-(ii+0.5*(1-i))*gridsp[0]);
        for(j=1;j>=-1;j-=2)
        {
            jjj+=j;
            dummy2 = dummy1 * j * (pt[1]-(jj+0.5*(1-j))*gridsp[1]);
            for(k=1;k>=-1;k-=2)
            {
                kkk+=k;
                factor = dummy2 * k * (pt[2]-(kk+0.5*(1-k))*gridsp[2]);
                msmul(stress,factor,gridres);
                m_add(sgrid[gmxLS_modulo(iii,nx)*nz*ny+gmxLS_modulo(jjj,ny)*nz+gmxLS_modulo(kkk,nz)],gridres,sgrid[gmxLS_modulo(iii,nx)*nz*ny+gmxLS_modulo(jjj,ny)*nz+gmxLS_modulo(kkk,nz)]);
            }
        }
    }
}

//----------------------------------------------------------------------------------------
//SPECIFIC DECOMPOSITIONS FOR 3,4 AND 5 PARTICLES


// Decompose 3-body potentials (angles)
void gmxLS_spread_n3(gmxLS_locals_grid_t * grid, rvec Fa, rvec Fb, rvec Fc, rvec Ra, rvec Rb, rvec Rc)
{

    //Counter
    int i;

    //************************************************************************************
    // UNIT vectors between particles
    rvec AB, AC, BC;

    // Distances
    real normAB,normAC,normBC;

    // (Covariant) Central Force decomposition
    real lab, lac, lbc;

    rvec Fij;
    //************************************************************************************

    //Dimension and number of particles
    int nDim = 3;
    int nPart = 3;

    //Number of rows and columns
    int nRow = nRow3, nCol = nCol3, nRHS = 1;

    //************************************************************************************
    // These are for the LAPACK dgelsd function
    real wkopt;
    real* work;
    int lwork = -1;
    int iwork[3*nRow3*0+11*nCol3];
    int info  =  0;
    int rank;
    real eps1 = epsLS;
    //************************************************************************************


    // Matrix of the system (12 equations x 6 unknowns)
    real M[nRow3*nCol3];
    // Vector, we want to solve M*x = b
    real b[nRow3], s[nCol3];

    //For MOP:
    rvec Rcenter1, Rcenter2;
    rvec Fcenter1, Fcenter2;
    int j;

    // If the force decomposition is cCFD or CFD
    if(grid->fdecomp == encCFD || grid->fdecomp == enCFD)
    {
        rvec_sub(Rb, Ra, AB);
        rvec_sub(Rc, Ra, AC);
        rvec_sub(Rc, Rb, BC);

        normAB=norm(AB);
        normAC=norm(AC);
        normBC=norm(BC);

        for ( i = 0; i < nCol*nRow; i++ )
        {
            M[i] = 0.0;
        }

        //Force on particle 1:
        M[nRow*0+0] = AB[0]; M[nRow*1+0] = AC[0];
        M[nRow*0+1] = AB[1]; M[nRow*1+1] = AC[1];
        M[nRow*0+2] = AB[2]; M[nRow*1+2] = AC[2];
        b[0] = Fa[0]; b[1] = Fa[1]; b[2] = Fa[2];

        //Force on particle 2:
        M[nRow*0+3] = -AB[0]; M[nRow*2+3] = BC[0];
        M[nRow*0+4] = -AB[1]; M[nRow*2+4] = BC[1];
        M[nRow*0+5] = -AB[2]; M[nRow*2+5] = BC[2];
        b[3] = Fb[0]; b[4] = Fb[1]; b[5] = Fb[2];

        //Force on particle 3:
        M[nRow*1+6] = -AC[0]; M[nRow*2+6] = -BC[0];
        M[nRow*1+7] = -AC[1]; M[nRow*2+7] = -BC[1];
        M[nRow*1+8] = -AC[2]; M[nRow*2+8] = -BC[2];
        b[6] = Fc[0]; b[7] = Fc[1]; b[8] = Fc[2];

        /* Query and allocate the optimal workspace */
        lwork = -1;
        F77_FUNC(dgelsd,DGELSD)(&nRow, &nCol, &nRHS, M, &nRow, b, &nRow, s, &eps1, &rank, &wkopt, &lwork, iwork, &info );
        lwork = (int)wkopt;
        work = (real*) malloc( lwork*sizeof(real) );
        /* Solve the equations A*X = B */
        // Least-Squares solution to the system 
        F77_FUNC(dgelsd,DGELSD)(&nRow, &nCol, &nRHS, M, &nRow, b, &nRow, s, &eps1, &rank, work, &lwork, iwork, &info );

        // Sum the 6 contributions to the stress

        lab = b[0];
        lac = b[1];
        lbc = b[2];

        Fij[0] = lab * AB[0]; Fij[1] = lab * AB[1]; Fij[2] = lab * AB[2];
        gmxLS_distribute_interaction(grid, Ra, Rb, Fij,0);
        Fij[0] = lac * AC[0]; Fij[1] = lac * AC[1]; Fij[2] = lac * AC[2];
        gmxLS_distribute_interaction(grid, Ra, Rc, Fij,0);
        Fij[0] = lbc * BC[0]; Fij[1] = lbc * BC[1]; Fij[2] = lbc * BC[2];
        gmxLS_distribute_interaction(grid, Rb, Rc, Fij,0);

        free(work);
    }
    else if (grid->fdecomp == enGLD)
    {
        Fij[0] = (Fa[XX]-Fb[XX])/3.0; Fij[1] = (Fa[YY]-Fb[YY])/3.0; Fij[2] = (Fa[ZZ]-Fb[ZZ])/3.0;
        gmxLS_distribute_interaction(grid, Ra, Rb, Fij,0);
        Fij[0] = (Fa[XX]-Fc[XX])/3.0; Fij[1] = (Fa[YY]-Fc[YY])/3.0; Fij[2] = (Fa[ZZ]-Fc[ZZ])/3.0;
        gmxLS_distribute_interaction(grid, Ra, Rc, Fij,0);
        Fij[0] = (Fb[XX]-Fc[XX])/3.0; Fij[1] = (Fb[YY]-Fc[YY])/3.0; Fij[2] = (Fb[ZZ]-Fc[ZZ])/3.0;
        gmxLS_distribute_interaction(grid, Rb, Rc, Fij,0);
    }
    else if (grid->fdecomp == enMOP)
    {
        double t;
        int skip;

        //All possible combinations:
        for ( i = 0; i < nDim; i++ )
        {
            // (12) - (3)
            // Center 1 (12)
            rvec_add(Ra, Rb, Rcenter1);
            svmul(0.5,Rcenter1,Rcenter1);

            //Center 2 (3)
            svmul(1.0,Rc,Rcenter2);
            rvec_add(Fa, Fb, Fcenter1);

            //Skip is 1 if the particles of the first group are further than the second geometric center in the i direction
            skip = 0;

            //Crossing with the plane in the i direction generated by the particle a
            t = (Ra[i]-Rcenter1[i])/(Rcenter2[i]-Rcenter1[i]);
            if (t > 1.0)
                skip = 1;
            else if ( t > 0 && skip == 0)
            {
                for ( j = 0; j < nDim; j++ )
                    Rcenter1[j] += t * (Rcenter2[j]-Rcenter1[j]);
            }

            //Crossing with the plane in the i direction generated by the particle b
            t = (Rb[i]-Rcenter1[i])/(Rcenter2[i]-Rcenter1[i]);
            if (t > 1.0)
                skip = 1;
            else if ( t > 0 && skip == 0)
            {
                for ( j = 0; j < nDim; j++ )
                    Rcenter1[j] += t * (Rcenter2[j]-Rcenter1[j]);
            }

            if (skip == 0)
                gmxLS_distribute_interaction(grid, Rcenter1, Rcenter2, Fcenter1,i+1);
        }

        for ( i = 0; i < nDim; i++ )
        {

            // (13) - (2)
            // Center 1 (12)
            rvec_add(Ra, Rc, Rcenter1);
            svmul(0.5,Rcenter1,Rcenter1);

            //Center 2 (2)
            svmul(1.0,Rb,Rcenter2);
            rvec_add(Fa, Fc, Fcenter1);

            //Skip is 1 if the particles of the first group are further than the second geometric center in the i direction
            skip = 0;

            //Crossing with the plane in the i direction generated by the particle a
            t = (Ra[i]-Rcenter1[i])/(Rcenter2[i]-Rcenter1[i]);
            if (t > 1.0)
                skip = 1;
            else if ( t > 0 && skip == 0)
            {
                for ( j = 0; j < nDim; j++ )
                    Rcenter1[j] += t * (Rcenter2[j]-Rcenter1[j]);
            }

            //Crossing with the plane in the i direction generated by the particle c
            t = (Rc[i]-Rcenter1[i])/(Rcenter2[i]-Rcenter1[i]);
            if (t > 1.0)
                skip = 1;
            else if ( t > 0 && skip == 0)
            {
                for ( j = 0; j < nDim; j++ )
                    Rcenter1[j] += t * (Rcenter2[j]-Rcenter1[j]);
            }

            if (skip == 0)
                gmxLS_distribute_interaction(grid, Rcenter1, Rcenter2, Fcenter1,i+1);
        }

        for ( i = 0; i < nDim; i++ )
        {

            // (23) - (1)
            // Center 1 (12)
            rvec_add(Rb, Rc, Rcenter1);
            svmul(0.5,Rcenter1,Rcenter1);

            //Center 2 (3)
            svmul(1.0,Ra,Rcenter2);
            rvec_add(Fb, Fc, Fcenter1);

            //Skip is 1 if the particles of the first group are further than the second geometric center in the i direction
            skip = 0;

            //Crossing with the plane in the i direction generated by the particle b
            t = (Rb[i]-Rcenter1[i])/(Rcenter2[i]-Rcenter1[i]);
            if (t > 1.0)
                skip = 1;
            else if ( t > 0 && skip == 0)
            {
                for ( j = 0; j < nDim; j++ )
                    Rcenter1[j] += t * (Rcenter2[j]-Rcenter1[j]);
            }

            //Crossing with the plane in the i direction generated by the particle c
            t = (Rc[i]-Rcenter1[i])/(Rcenter2[i]-Rcenter1[i]);
            if (t > 1.0)
                skip = 1;
            else if ( t > 0 && skip == 0)
            {
                for ( j = 0; j < nDim; j++ )
                    Rcenter1[j] += t * (Rcenter2[j]-Rcenter1[j]);
            }

            if (skip==0)
                gmxLS_distribute_interaction(grid, Rcenter1, Rcenter2, Fcenter1,i+1);
        }

    }
}

// Hybrid Decomposition for angle potentials
void gmxLS_spread_n3_HD(gmxLS_locals_grid_t * grid, rvec F1, rvec F2, rvec F3, rvec r1, rvec r2, rvec r3)
{
  int i;

  rvec dr12, dr23, dr13;
  rvec dr41, dr42, dr43;
  rvec r4;

  rvec Fij;

  // Matrix of the system (12 equations x 6 unknowns)
  real M[nRow3_HD * nCol3_HD] = {0.0};
  // Vector, we want to solve M*x = b
  real b[nRow3_HD] = {0.0}, s[nCol3_HD] = {0.0};

  int nRow = nRow3_HD, nCol = nCol3_HD, nRHS = 1;

  real wkopt;
  real* work;
  int lwork = -1;
  int iwork[3 * nRow3_HD * 0 + 11 * nCol3_HD];
  int info  =  0;
  int rank;
  real eps1 = epsLS;

  rvec_sub(r2, r1, dr12); // r2 - r1
  rvec_sub(r3, r2, dr23); // r3 - r2
  rvec_sub(r3, r1, dr13); // r3 - r1
\
  // NOTE: In order to verify the original code of CFD,
  //       we also implement CFD.
  if (grid->fdecomp == encCFD || grid->fdecomp == enCFD) {
    rvec_add(r1, r3, r4);
    svmul(0.5, r4, r4);
  } else if (grid->fdecomp == enFCD) {
    gmxLS_get_OFC(r1, r2, dr12, F2, r4);
  } else if (grid->fdecomp == enHD_GM) {
    gmxLS_get_SDM_gmin(r1, r2, r3, dr12, dr23, r4);
  } else if (grid->fdecomp == enHD_LM) {
    gmxLS_get_SDM_lmin(r1, r2, r3, dr12, dr23, r4);
  } else {
    fprintf(stderr, "Unknown force decomposition mode at %s %d\n", __FILE__, __LINE__);
    exit(-1);
  }

  rvec_sub(r1, r4, dr41);
  rvec_sub(r2, r4, dr42);
  rvec_sub(r3, r4, dr43);

  // lapack solver
  M[nRow3_HD * 0 + 0] = -dr12[0]; M[nRow3_HD * 1 + 0] = dr41[0];
  M[nRow3_HD * 0 + 1] = -dr12[1]; M[nRow3_HD * 1 + 1] = dr41[1];
  M[nRow3_HD * 0 + 2] = -dr12[2]; M[nRow3_HD * 1 + 2] = dr41[2];
  b[0] = F1[0]; b[1] = F1[1]; b[2] = F1[2];

  M[nRow3_HD * 0 + 3] =  dr12[0]; M[nRow3_HD * 2 + 3] = -dr23[0]; M[nRow3_HD * 3 + 3] = dr42[0];
  M[nRow3_HD * 0 + 4] =  dr12[1]; M[nRow3_HD * 2 + 4] = -dr23[1]; M[nRow3_HD * 3 + 4] = dr42[1];
  M[nRow3_HD * 0 + 5] =  dr12[2]; M[nRow3_HD * 2 + 5] = -dr23[2]; M[nRow3_HD * 3 + 5] = dr42[2];
  b[3] = F2[0]; b[4] = F2[1]; b[5] = F2[2];

  M[nRow3_HD * 2 + 6] = dr23[0]; M[nRow3_HD * 4 + 6] = dr43[0];
  M[nRow3_HD * 2 + 7] = dr23[1]; M[nRow3_HD * 4 + 7] = dr43[1];
  M[nRow3_HD * 2 + 8] = dr23[2]; M[nRow3_HD * 4 + 8] = dr43[2];
  b[6] = F3[0]; b[7] = F3[1]; b[8] = F3[2];

  M[nRow3_HD * 1 + 9]  = -dr41[0]; M[nRow3_HD * 3 + 9]  = -dr42[0]; M[nRow3_HD * 4 + 9]  = -dr43[0];
  M[nRow3_HD * 1 + 10] = -dr41[1]; M[nRow3_HD * 3 + 10] = -dr42[1]; M[nRow3_HD * 4 + 10] = -dr43[1];
  M[nRow3_HD * 1 + 11] = -dr41[2]; M[nRow3_HD * 3 + 11] = -dr42[2]; M[nRow3_HD * 4 + 11] = -dr43[2];
  b[9] = 0.0;  b[10] = 0.0; b[11] = 0.0;

  lwork = -1;
  F77_FUNC(dgelsd, DGELSD)(&nRow, &nCol, &nRHS, M, &nRow, b, &nRow, s, &eps1, &rank, &wkopt, &lwork, iwork, &info);
  lwork = (int) wkopt;
  work = (real*) malloc(lwork * sizeof(real));

  F77_FUNC(dgelsd, DGELSD)(&nRow, &nCol, &nRHS, M, &nRow, b, &nRow, s, &eps1, &rank, work, &lwork, iwork, &info);

  // spread stress
  Fij[0] = b[0] * (-dr12[0]); Fij[1] = b[0] * (-dr12[1]); Fij[2] = b[0] * (-dr12[2]);
  gmxLS_distribute_interaction(grid, r2, r1, Fij, 0);
  Fij[0] = b[1] * dr41[0]; Fij[1] = b[1] * dr41[1]; Fij[2] = b[1] * dr41[2];
  gmxLS_distribute_interaction(grid, r4, r1, Fij, 0);
  Fij[0] = b[2] * (-dr23[0]); Fij[1] = b[2] * (-dr23[1]); Fij[2] = b[2] * (-dr23[2]);
  gmxLS_distribute_interaction(grid, r3, r2, Fij, 0);
  Fij[0] = b[3] * dr42[0]; Fij[1] = b[3] * dr42[1]; Fij[2] = b[3] * dr42[2];
  gmxLS_distribute_interaction(grid, r4, r2, Fij, 0);
  Fij[0] = b[4] * dr43[0]; Fij[1] = b[4] * dr43[1]; Fij[2] = b[4] * dr43[2];
  gmxLS_distribute_interaction(grid, r4, r3, Fij, 0);

  free(work);
}

void gmxLS_get_SDM_min(const rvec r1, const rvec r2, const rvec r3, const rvec dr12, const rvec dr23, const real sign, rvec ret)
{
  real dr12_norm, dr23_norm, dr24_norm, cos123;
  rvec dr21_hat, dr23_hat, dr24_hat;
  rvec dr24;

  dr12_norm = norm(dr12);
  dr23_norm = norm(dr23);
  dr24_norm = sqrt(dr12_norm * dr23_norm);

  cos123 = -1.0 * iprod(dr12, dr23) / (dr12_norm * dr23_norm);
  assert(dr12_norm >= (dr23_norm * (1.0 + cos123) * 0.5));
  assert(dr23_norm >= (dr12_norm * (1.0 + cos123) * 0.5));

  svmul(-1.0 / dr12_norm, dr12, dr21_hat);
  svmul( 1.0 / dr23_norm, dr23, dr23_hat);

  rvec_add(dr21_hat, dr23_hat, dr24_hat);
  unitv(dr24_hat, dr24_hat);

  svmul(dr24_norm * sign, dr24_hat, dr24);
  rvec_add(r2, dr24, ret);
}

void gmxLS_get_SDM_gmin(const rvec r1, const rvec r2, const rvec r3, const rvec dr12, const rvec dr23, rvec ret)
{
  gmxLS_get_SDM_min(r1, r2, r3, dr12, dr23, 1.0, ret);
}

void gmxLS_get_SDM_lmin(const rvec r1, const rvec r2, const rvec r3, const rvec dr12, const rvec dr23, rvec ret)
{
  gmxLS_get_SDM_min(r1, r2, r3, dr12, dr23, -1.0, ret);
}

void gmxLS_get_OFC(const rvec r1, const rvec r2, const rvec dr12, const rvec F2, rvec ret)
{
  rvec dr21;
  rvec offset;
  real dr21_norm2, F2_dr21;
  /* r2 + F2 * (dr21 * dr21 / (F2 * dr21)) */

  svmul(-1.0, dr12, dr21);

  dr21_norm2 = norm2(dr21);
  F2_dr21 = iprod(F2, dr21);
  svmul(dr21_norm2 / F2_dr21, F2, offset);
  rvec_add(r2, offset, ret);
}

// Decompose Settle
void gmxLS_spread_settle(gmxLS_locals_grid_t * grid, rvec Fa, rvec Fb, rvec Fc, rvec Ra, rvec Rb, rvec Rc)
{
    //Counter
    int i;

    //************************************************************************************
    // UNIT vectors between particles
    rvec AB, AC, BC;

    // Distances
    real normAB,normAC,normBC;

    // (Covariant) Central Force decomposition
    real lab, lac, lbc;

    rvec Fij;
    //************************************************************************************

    //Dimension and number of particles
    int nDim = 3;
    int nPart = 3;

    //Number of rows and columns
    int nRow = nRow3, nCol = nCol3, nRHS = 1;

    //************************************************************************************
    // These are for the LAPACK dgelsd function
    real wkopt;
    real* work;
    int lwork = -1;
    int iwork[3*nRow3*0+11*nCol3];
    int info  =  0;
    int rank;
    real eps1 = epsLS;
    //************************************************************************************

    // Matrix of the system (12 equations x 6 unknowns)
    real M[nRow3*nCol3];
    // Vector, we want to solve M*x = b
    real b[nRow3], s[nCol3];

    if (grid->fdecomp == encCFD || grid->fdecomp == enCFD || grid->fdecomp == enGLD || grid->fdecomp == enMOP)
    {
        rvec_sub(Rb, Ra, AB);
        rvec_sub(Rc, Ra, AC);
        rvec_sub(Rc, Rb, BC);

        normAB=norm(AB);
        normAC=norm(AC);
        normBC=norm(BC);

        for ( i = 0; i < nCol*nRow; i++ )
        {
            M[i] = 0.0;
        }

        //Force on particle 1:
        M[nRow*0+0] = AB[0]; M[nRow*1+0] = AC[0];
        M[nRow*0+1] = AB[1]; M[nRow*1+1] = AC[1];
        M[nRow*0+2] = AB[2]; M[nRow*1+2] = AC[2];
        b[0] = Fa[0]; b[1] = Fa[1]; b[2] = Fa[2];

        //Force on particle 2:
        M[nRow*0+3] = -AB[0]; M[nRow*2+3] = BC[0];
        M[nRow*0+4] = -AB[1]; M[nRow*2+4] = BC[1];
        M[nRow*0+5] = -AB[2]; M[nRow*2+5] = BC[2];
        b[3] = Fb[0]; b[4] = Fb[1]; b[5] = Fb[2];

        //Force on particle 3:
        M[nRow*1+6] = -AC[0]; M[nRow*2+6] = -BC[0];
        M[nRow*1+7] = -AC[1]; M[nRow*2+7] = -BC[1];
        M[nRow*1+8] = -AC[2]; M[nRow*2+8] = -BC[2];
        b[6] = Fc[0]; b[7] = Fc[1]; b[8] = Fc[2];

        /* Query and allocate the optimal workspace */
        lwork = -1;
        F77_FUNC(dgelsd,DGELSD)(&nRow, &nCol, &nRHS, M, &nRow, b, &nRow, s, &eps1, &rank, &wkopt, &lwork, iwork, &info );
        lwork = (int)wkopt;
        work = (real*) malloc( lwork*sizeof(real) );
        /* Solve the equations A*X = B */
        // Least-Squares solution to the system
        F77_FUNC(dgelsd,DGELSD)(&nRow, &nCol, &nRHS, M, &nRow, b, &nRow, s, &eps1, &rank, work, &lwork, iwork, &info );

        // Sum the 3 contributions to the stress
        lab = b[0];
        lac = b[1];
        lbc = b[2];

        Fij[0] = lab * AB[0]; Fij[1] = lab * AB[1]; Fij[2] = lab * AB[2];
        gmxLS_distribute_interaction(grid, Ra, Rb, Fij,0);
        Fij[0] = lac * AC[0]; Fij[1] = lac * AC[1]; Fij[2] = lac * AC[2];
        gmxLS_distribute_interaction(grid, Ra, Rc, Fij,0);
        Fij[0] = lbc * BC[0]; Fij[1] = lbc * BC[1]; Fij[2] = lbc * BC[2];
        gmxLS_distribute_interaction(grid, Rb, Rc, Fij,0);

        free(work);
    }
}

// Decompose 4-body potentials (dihedrals)
void gmxLS_spread_n4(gmxLS_locals_grid_t * grid, rvec Fa, rvec Fb, rvec Fc, rvec Fd, rvec Ra, rvec Rb, rvec Rc, rvec Rd)
{
    //Counter
    int i;

    //************************************************************************************
    // UNIT vectors between particles
    rvec AB, AC, AD, BC, BD, CD;

    // Distances
    real normAB,normAC,normAD,normBC,normBD,normCD;

    // (Covariant) Central Force decomposition
    real lab, lac, lad, lbc, lbd, lcd;

    rvec Fij;
    //************************************************************************************

    //Dimension and number of particles
    int nDim = 3;
    int nPart = 4;

    //Number of rows and columns
    int nRow = nRow4, nCol = nCol4, nRHS = 1;

    //************************************************************************************
    // These are for the LAPACK dgelsd function
    real wkopt;
    real* work;
    int lwork = -1;
    int iwork[3*nRow4*0+11*nCol4];
    int info  =  0;
    int rank;
    real eps1 = epsLS;
    //************************************************************************************

    // Matrix of the system (12 equations x 6 unknowns)
    real M[nRow4*nCol4];
    // Vector, we want to solve M*x = b
    real b[nRow4], s[nCol4];

    //For MOP:
    rvec Rcenter1, Rcenter2;
    rvec Fcenter1, Fcenter2;
    int j;

    // If the force decomposition is cCFD or CFD
    if(grid->fdecomp == encCFD || grid->fdecomp == enCFD)
    {
        rvec_sub(Rb, Ra, AB);
        rvec_sub(Rc, Ra, AC);
        rvec_sub(Rd, Ra, AD);
        rvec_sub(Rc, Rb, BC);
        rvec_sub(Rd, Rb, BD);
        rvec_sub(Rd, Rc, CD);

        normAB=norm(AB);
        normAC=norm(AC);
        normAD=norm(AD);
        normBC=norm(BC);
        normBD=norm(BD);
        normCD=norm(CD);

        for ( i = 0; i < nCol*nRow; i++ )
        {
            M[i] = 0.0;
        }
        //Force on particle 1:
        M[nRow*0+0] = AB[0]; M[nRow*1+0] = AC[0]; M[nRow*2+0] = AD[0];
        M[nRow*0+1] = AB[1]; M[nRow*1+1] = AC[1]; M[nRow*2+1] = AD[1];
        M[nRow*0+2] = AB[2]; M[nRow*1+2] = AC[2]; M[nRow*2+2] = AD[2];
        b[0] = Fa[0]; b[1] = Fa[1]; b[2] = Fa[2];

        //Force on particle 2:
        M[nRow*0+3] = -AB[0]; M[nRow*3+3] = BC[0]; M[nRow*4+3] = BD[0];
        M[nRow*0+4] = -AB[1]; M[nRow*3+4] = BC[1]; M[nRow*4+4] = BD[1];
        M[nRow*0+5] = -AB[2]; M[nRow*3+5] = BC[2]; M[nRow*4+5] = BD[2];
        b[3] = Fb[0]; b[4] = Fb[1]; b[5] = Fb[2];

        //Force on particle 3:
        M[nRow*1+6] = -AC[0]; M[nRow*3+6] = -BC[0]; M[nRow*5+6] = CD[0];
        M[nRow*1+7] = -AC[1]; M[nRow*3+7] = -BC[1]; M[nRow*5+7] = CD[1];
        M[nRow*1+8] = -AC[2]; M[nRow*3+8] = -BC[2]; M[nRow*5+8] = CD[2];
        b[6] = Fc[0]; b[7] = Fc[1]; b[8] = Fc[2];

        //Force on particle 4:
        M[nRow*2+9]  = -AD[0]; M[nRow*4+9] = -BD[0]; M[nRow*5+9] = -CD[0];
        M[nRow*2+10] = -AD[1]; M[nRow*4+10] = -BD[1]; M[nRow*5+10] = -CD[1];
        M[nRow*2+11] = -AD[2]; M[nRow*4+11] = -BD[2]; M[nRow*5+11] = -CD[2];
        b[9] = Fd[0]; b[10] = Fd[1]; b[11] = Fd[2];

        /* Query and allocate the optimal workspace */
        lwork = -1;
        F77_FUNC(dgelsd,DGELSD)(&nRow, &nCol, &nRHS, M, &nRow, b, &nRow, s, &eps1, &rank, &wkopt, &lwork, iwork, &info );
        lwork = (int)wkopt;
        work = (real*) malloc( lwork*sizeof(real) );
        /* Solve the equations A*X = B */
        // Least-Squares solution to the system
        F77_FUNC(dgelsd,DGELSD)(&nRow, &nCol, &nRHS, M, &nRow, b, &nRow, s, &eps1, &rank, work, &lwork, iwork, &info );

        // Sum the 6 contributions to the stress

        lab = b[0];
        lac = b[1];
        lad = b[2];
        lbc = b[3];
        lbd = b[4];
        lcd = b[5];

        Fij[0] = lab * AB[0]; Fij[1] = lab * AB[1]; Fij[2] = lab * AB[2];
        gmxLS_distribute_interaction(grid, Ra, Rb, Fij,0);
        Fij[0] = lac * AC[0]; Fij[1] = lac * AC[1]; Fij[2] = lac * AC[2];
        gmxLS_distribute_interaction(grid, Ra, Rc, Fij,0);
        Fij[0] = lad * AD[0]; Fij[1] = lad * AD[1]; Fij[2] = lad * AD[2];
        gmxLS_distribute_interaction(grid, Ra, Rd, Fij,0);
        Fij[0] = lbc * BC[0]; Fij[1] = lbc * BC[1]; Fij[2] = lbc * BC[2];
        gmxLS_distribute_interaction(grid, Rb, Rc, Fij,0);
        Fij[0] = lbd * BD[0]; Fij[1] = lbd * BD[1]; Fij[2] = lbd * BD[2];
        gmxLS_distribute_interaction(grid, Rb, Rd, Fij,0);
        Fij[0] = lcd * CD[0]; Fij[1] = lcd * CD[1]; Fij[2] = lcd * CD[2];
        gmxLS_distribute_interaction(grid, Rc, Rd, Fij,0);

        free(work);
    }
    else if (grid->fdecomp == enGLD)
    {
        Fij[0] = (Fa[XX]-Fb[XX])/4.0; Fij[1] = (Fa[YY]-Fb[YY])/4.0; Fij[2] = (Fa[ZZ]-Fb[ZZ])/4.0;
        gmxLS_distribute_interaction(grid, Ra, Rb, Fij,0);
        Fij[0] = (Fa[XX]-Fc[XX])/4.0; Fij[1] = (Fa[YY]-Fc[YY])/4.0; Fij[2] = (Fa[ZZ]-Fc[ZZ])/4.0;
        gmxLS_distribute_interaction(grid, Ra, Rc, Fij,0);
        Fij[0] = (Fa[XX]-Fd[XX])/4.0; Fij[1] = (Fa[YY]-Fd[YY])/4.0; Fij[2] = (Fa[ZZ]-Fd[ZZ])/4.0;
        gmxLS_distribute_interaction(grid, Ra, Rd, Fij,0);
        Fij[0] = (Fb[XX]-Fc[XX])/4.0; Fij[1] = (Fb[YY]-Fc[YY])/4.0; Fij[2] = (Fb[ZZ]-Fc[ZZ])/4.0;
        gmxLS_distribute_interaction(grid, Rb, Rc, Fij,0);
        Fij[0] = (Fb[XX]-Fd[XX])/4.0; Fij[1] = (Fb[YY]-Fd[YY])/4.0; Fij[2] = (Fb[ZZ]-Fd[ZZ])/4.0;
        gmxLS_distribute_interaction(grid, Rb, Rd, Fij,0);
        Fij[0] = (Fc[XX]-Fd[XX])/4.0; Fij[1] = (Fc[YY]-Fd[YY])/4.0; Fij[2] = (Fc[ZZ]-Fd[ZZ])/4.0;
        gmxLS_distribute_interaction(grid, Rc, Rd, Fij,0);
    }
    else if (grid->fdecomp == enMOP)
    {
        double t;
        int skip;

        //All possible combinations:
        for ( i = 0; i < nDim; i++ )
        {
            // (123) - (4)
            // Center 1 (123)
            rvec_add(Ra, Rb, Rcenter1);
            rvec_add(Rcenter1, Rc, Rcenter1);
            svmul(0.333333333333,Rcenter1,Rcenter1);

            //Center 2 (4)
            svmul(1.0,Rd,Rcenter2);

            rvec_add(Fa, Fb, Fcenter1);
            rvec_add(Fcenter1, Fc, Fcenter1);

            //Skip is 1 if the particles of the first group are further than the second geometric center in the i direction
            skip = 0;

            //Crossing with the plane in the i direction generated by the particle a
            t = (Ra[i]-Rcenter1[i])/(Rcenter2[i]-Rcenter1[i]);
            if (t > 1.0)
                skip = 1;
            else if ( t > 0 && skip == 0)
            {
                for ( j = 0; j < nDim; j++ )
                    Rcenter1[j] += t * (Rcenter2[j]-Rcenter1[j]);
            }

            //Crossing with the plane in the i direction generated by the particle b
            t = (Rb[i]-Rcenter1[i])/(Rcenter2[i]-Rcenter1[i]);
            if (t > 1.0)
                skip = 1;
            else if ( t > 0 && skip == 0)
            {
                for ( j = 0; j < nDim; j++ )
                    Rcenter1[j] += t * (Rcenter2[j]-Rcenter1[j]);
            }

            //Crossing with the plane in the i direction generated by the particle b
            t = (Rc[i]-Rcenter1[i])/(Rcenter2[i]-Rcenter1[i]);
            if (t > 1.0)
                skip = 1;
            else if ( t > 0 && skip == 0)
            {
                for ( j = 0; j < nDim; j++ )
                    Rcenter1[j] += t * (Rcenter2[j]-Rcenter1[j]);
            }

            if (skip == 0)
                gmxLS_distribute_interaction(grid, Rcenter1, Rcenter2, Fcenter1,i+1);
        }

        for ( i = 0; i < nDim; i++ )
        {

            // (124) - (3)
            // Center 1 (124)
            rvec_add(Ra, Rb, Rcenter1);
            rvec_add(Rcenter1, Rd, Rcenter1);
            svmul(0.333333333333,Rcenter1,Rcenter1);

            //Center 2 (3)
            svmul(1.0,Rc,Rcenter2);

            rvec_add(Fa, Fb, Fcenter1);
            rvec_add(Fcenter1, Fd, Fcenter1);

            //Skip is 1 if the particles of the first group are further than the second geometric center in the i direction
            skip = 0;

            //Crossing with the plane in the i direction generated by the particle a
            t = (Ra[i]-Rcenter1[i])/(Rcenter2[i]-Rcenter1[i]);
            if (t > 1.0)
                skip = 1;
            else if ( t > 0 && skip == 0)
            {
                for ( j = 0; j < nDim; j++ )
                    Rcenter1[j] += t * (Rcenter2[j]-Rcenter1[j]);
            }

            //Crossing with the plane in the i direction generated by the particle b
            t = (Rb[i]-Rcenter1[i])/(Rcenter2[i]-Rcenter1[i]);
            if (t > 1.0)
                skip = 1;
            else if ( t > 0 && skip == 0)
            {
                for ( j = 0; j < nDim; j++ )
                    Rcenter1[j] += t * (Rcenter2[j]-Rcenter1[j]);
            }

            //Crossing with the plane in the i direction generated by the particle b
            t = (Rd[i]-Rcenter1[i])/(Rcenter2[i]-Rcenter1[i]);
            if (t > 1.0)
                skip = 1;
            else if ( t > 0 && skip == 0)
            {
                for ( j = 0; j < nDim; j++ )
                    Rcenter1[j] += t * (Rcenter2[j]-Rcenter1[j]);
            }

            if (skip == 0)
                gmxLS_distribute_interaction(grid, Rcenter1, Rcenter2, Fcenter1,i+1);
        }

        for ( i = 0; i < nDim; i++ )
        {

            // (134) - (2)
            // Center 1 (134)
            rvec_add(Ra, Rc, Rcenter1);
            rvec_add(Rcenter1, Rd, Rcenter1);
            svmul(0.333333333333,Rcenter1,Rcenter1);

            //Center 2 (2)
            svmul(1.0,Rb,Rcenter2);

            rvec_add(Fa, Fc, Fcenter1);
            rvec_add(Fcenter1, Fd, Fcenter1);

            //Skip is 1 if the particles of the first group are further than the second geometric center in the i direction
            skip = 0;

            //Crossing with the plane in the i direction generated by the particle a
            t = (Ra[i]-Rcenter1[i])/(Rcenter2[i]-Rcenter1[i]);
            if (t > 1.0)
                skip = 1;
            else if ( t > 0 && skip == 0)
            {
                for ( j = 0; j < nDim; j++ )
                    Rcenter1[j] += t * (Rcenter2[j]-Rcenter1[j]);
            }

            //Crossing with the plane in the i direction generated by the particle b
            t = (Rc[i]-Rcenter1[i])/(Rcenter2[i]-Rcenter1[i]);
            if (t > 1.0)
                skip = 1;
            else if ( t > 0 && skip == 0)
            {
                for ( j = 0; j < nDim; j++ )
                    Rcenter1[j] += t * (Rcenter2[j]-Rcenter1[j]);
            }

            //Crossing with the plane in the i direction generated by the particle b
            t = (Rd[i]-Rcenter1[i])/(Rcenter2[i]-Rcenter1[i]);
            if (t > 1.0)
                skip = 1;
            else if ( t > 0 && skip == 0)
            {
                for ( j = 0; j < nDim; j++ )
                    Rcenter1[j] += t * (Rcenter2[j]-Rcenter1[j]);
            }

            if (skip == 0)
                gmxLS_distribute_interaction(grid, Rcenter1, Rcenter2, Fcenter1,i+1);
        }

        for ( i = 0; i < nDim; i++ )
        {

            // (234) - (1)
            // Center 1 (234)
            rvec_add(Rb, Rc, Rcenter1);
            rvec_add(Rcenter1, Rd, Rcenter1);
            svmul(0.333333333333,Rcenter1,Rcenter1);

            //Center 2 (1)
            svmul(1.0,Ra,Rcenter2);

            rvec_add(Fb, Fc, Fcenter1);
            rvec_add(Fcenter1, Fd, Fcenter1);

            //Skip is 1 if the particles of the first group are further than the second geometric center in the i direction
            skip = 0;

            //Crossing with the plane in the i direction generated by the particle a
            t = (Rb[i]-Rcenter1[i])/(Rcenter2[i]-Rcenter1[i]);
            if (t > 1.0)
                skip = 1;
            else if ( t > 0 && skip == 0)
            {
                for ( j = 0; j < nDim; j++ )
                    Rcenter1[j] += t * (Rcenter2[j]-Rcenter1[j]);
            }

            //Crossing with the plane in the i direction generated by the particle b
            t = (Rc[i]-Rcenter1[i])/(Rcenter2[i]-Rcenter1[i]);
            if (t > 1.0)
                skip = 1;
            else if ( t > 0 && skip == 0)
            {
                for ( j = 0; j < nDim; j++ )
                    Rcenter1[j] += t * (Rcenter2[j]-Rcenter1[j]);
            }

            //Crossing with the plane in the i direction generated by the particle b
            t = (Rd[i]-Rcenter1[i])/(Rcenter2[i]-Rcenter1[i]);
            if (t > 1.0)
                skip = 1;
            else if ( t > 0 && skip == 0)
            {
                for ( j = 0; j < nDim; j++ )
                    Rcenter1[j] += t * (Rcenter2[j]-Rcenter1[j]);
            }

            if (skip == 0)
                gmxLS_distribute_interaction(grid, Rcenter1, Rcenter2, Fcenter1,i+1);
        }

        //All possible combinations:
        for ( i = 0; i < nDim; i++ )
        {

            // (12) - (34)
            // Center 1 (12)
            rvec_add(Ra, Rb, Rcenter1);
            svmul(0.5,Rcenter1,Rcenter1);

            // Center 2 (34)
            rvec_add(Rc, Rd, Rcenter2);
            svmul(0.5,Rcenter2,Rcenter2);

            rvec_add(Fa, Fb, Fcenter1);

            //Skip is 1 if the particles of the first group are further than the second geometric center in the i direction
            skip = 0;

            //Crossing with the plane in the i direction generated by the particle a
            t = (Ra[i]-Rcenter1[i])/(Rcenter2[i]-Rcenter1[i]);
            if (t > 1.0)
                skip = 1;
            else if ( t > 0 && skip == 0)
            {
                for ( j = 0; j < nDim; j++ )
                    Rcenter1[j] += t * (Rcenter2[j]-Rcenter1[j]);
            }

            //Crossing with the plane in the i direction generated by the particle b
            t = (Rb[i]-Rcenter1[i])/(Rcenter2[i]-Rcenter1[i]);
            if (t > 1.0)
                skip = 1;
            else if ( t > 0 && skip == 0)
            {
                for ( j = 0; j < nDim; j++ )
                    Rcenter1[j] += t * (Rcenter2[j]-Rcenter1[j]);
            }

            //Crossing with the plane in the i direction generated by the particle c
            t = (Rc[i]-Rcenter1[i])/(Rcenter2[i]-Rcenter1[i]);
            if ( t < 0.0 )
                skip = 1;
            else if ( t < 1.0 && skip == 0)
            {
                for ( j = 0; j < nDim; j++ )
                    Rcenter2[j] += (1.0-t) * (Rcenter1[j]-Rcenter2[j]);
            }

            //Crossing with the plane in the i direction generated by the particle d
            t = (Rd[i]-Rcenter1[i])/(Rcenter2[i]-Rcenter1[i]);
            if ( t < 0.0 )
                skip = 1;
            else if ( t < 1.0 && skip == 0)
            {
                for ( j = 0; j < nDim; j++ )
                    Rcenter2[j] += (1.0-t) * (Rcenter1[j]-Rcenter2[j]);
            }

            if (skip == 0)
                gmxLS_distribute_interaction(grid, Rcenter1, Rcenter2, Fcenter1,i+1);
        }
        for ( i = 0; i < nDim; i++ )
        {
            // (13) - (24)
            // Center 1 (13)
            rvec_add(Ra, Rc, Rcenter1);
            svmul(0.5,Rcenter1,Rcenter1);

            // Center 2 (34)
            rvec_add(Rb, Rd, Rcenter2);
            svmul(0.5,Rcenter2,Rcenter2);

            rvec_add(Fa, Fc, Fcenter1);

            //Skip is 1 if the particles of the first group are further than the second geometric center in the i direction
            skip = 0;

            //Crossing with the plane in the i direction generated by the particle a
            t = (Ra[i]-Rcenter1[i])/(Rcenter2[i]-Rcenter1[i]);
            if (t > 1.0)
                skip = 1;
            else if ( t > 0 && skip == 0)
            {
                for ( j = 0; j < nDim; j++ )
                    Rcenter1[j] += t * (Rcenter2[j]-Rcenter1[j]);
            }

            //Crossing with the plane in the i direction generated by the particle b
            t = (Rc[i]-Rcenter1[i])/(Rcenter2[i]-Rcenter1[i]);
            if (t > 1.0)
                skip = 1;
            else if ( t > 0 && skip == 0)
            {
                for ( j = 0; j < nDim; j++ )
                    Rcenter1[j] += t * (Rcenter2[j]-Rcenter1[j]);
            }

            //Crossing with the plane in the i direction generated by the particle c
            t = (Rb[i]-Rcenter1[i])/(Rcenter2[i]-Rcenter1[i]);
            if ( t < 0.0 )
                skip = 1;
            else if ( t < 1.0 && skip == 0)
            {
                for ( j = 0; j < nDim; j++ )
                    Rcenter2[j] += (1.0-t) * (Rcenter1[j]-Rcenter2[j]);
            }

            //Crossing with the plane in the i direction generated by the particle d
            t = (Rd[i]-Rcenter1[i])/(Rcenter2[i]-Rcenter1[i]);
            if ( t < 0.0 )
                skip = 1;
            else if ( t < 1.0 && skip == 0)
            {
                for ( j = 0; j < nDim; j++ )
                    Rcenter2[j] += (1.0-t) * (Rcenter1[j]-Rcenter2[j]);
            }

            if (skip == 0)
                gmxLS_distribute_interaction(grid, Rcenter1, Rcenter2, Fcenter1,i+1);
        }
        for ( i = 0; i < nDim; i++ )
        {
            // (14) - (23)
            // Center 1 (14)
            rvec_add(Ra, Rd, Rcenter1);
            svmul(0.5,Rcenter1,Rcenter1);

            // Center 2 (34)
            rvec_add(Rb, Rc, Rcenter2);
            svmul(0.5,Rcenter2,Rcenter2);

            rvec_add(Fa, Fd, Fcenter1);

            //Skip is 1 if the particles of the first group are further than the second geometric center in the i direction
            skip = 0;

            //Crossing with the plane in the i direction generated by the particle a
            t = (Ra[i]-Rcenter1[i])/(Rcenter2[i]-Rcenter1[i]);
            if (t > 1.0)
                skip = 1;
            else if ( t > 0 && skip == 0)
            {
                for ( j = 0; j < nDim; j++ )
                    Rcenter1[j] += t * (Rcenter2[j]-Rcenter1[j]);
            }

            //Crossing with the plane in the i direction generated by the particle b
            t = (Rd[i]-Rcenter1[i])/(Rcenter2[i]-Rcenter1[i]);
            if (t > 1.0)
                skip = 1;
            else if ( t > 0 && skip == 0)
            {
                for ( j = 0; j < nDim; j++ )
                    Rcenter1[j] += t * (Rcenter2[j]-Rcenter1[j]);
            }

            //Crossing with the plane in the i direction generated by the particle c
            t = (Rb[i]-Rcenter1[i])/(Rcenter2[i]-Rcenter1[i]);
            if ( t < 0.0 )
                skip = 1;
            else if ( t < 1.0 && skip == 0)
            {
                for ( j = 0; j < nDim; j++ )
                    Rcenter2[j] += (1.0-t) * (Rcenter1[j]-Rcenter2[j]);
            }

            //Crossing with the plane in the i direction generated by the particle d
            t = (Rc[i]-Rcenter1[i])/(Rcenter2[i]-Rcenter1[i]);
            if ( t < 0.0 )
                skip = 1;
            else if ( t < 1.0 && skip == 0)
            {
                for ( j = 0; j < nDim; j++ )
                    Rcenter2[j] += (1.0-t) * (Rcenter1[j]-Rcenter2[j]);
            }

            if (skip == 0)
                gmxLS_distribute_interaction(grid, Rcenter1, Rcenter2, Fcenter1,i+1);
        }
    }
}


// Decompose 5-body potentials (CMAP)
void gmxLS_spread_n5(gmxLS_locals_grid_t * grid, rvec Fa, rvec Fb, rvec Fc, rvec Fd, rvec Fe, rvec Ra, rvec Rb, rvec Rc, rvec Rd, rvec Re)
{
    //Counter
    int i;

    //************************************************************************************
    // UNIT vectors between particles
    rvec AB, AC, AD, AE, BC, BD, BE, CD, CE, DE;

    // Distances
    real normAB,normAC,normAD,normAE,normBC,normBD,normBE,normCD,normCE,normDE;

    // (Covariant) Central Force decomposition
    real lab, lac, lad, lae, lbc, lbd, lbe, lcd, lce, lde;

    rvec Fij;
    //************************************************************************************

    //Dimension and number of particles
    int nDim = 3;
    int nPart = 5;

    //Number of rows and columns
    int nRow = nRow5, nCol = nCol5, nRHS = 1;

    //************************************************************************************
    // These are for the LAPACK dgelsd function
    real wkopt;
    real* work;
    int lwork = -1;
    int iwork[3*nRow5*0+11*nCol5];
    int info  =  0;
    int rank;
    real eps1 = epsLS;
    //************************************************************************************

    // Matrix of the system (15 equations x 10 unknowns)
    real M[nRow5*nCol5];
    // Vector, we want to solve M*x = b
    real b[nRow5], s[nCol5], CaleyMengerNormal[nCol5];
    // Scalar product of the Normal and the initial CFD
    real prod;

    // If the force decomposition is cCFD or CFD
    if(grid->fdecomp == encCFD || grid->fdecomp == enCFD)
    {
        rvec_sub(Rb, Ra, AB);
        rvec_sub(Rc, Ra, AC);
        rvec_sub(Rd, Ra, AD);
        rvec_sub(Re, Ra, AE);
        rvec_sub(Rc, Rb, BC);
        rvec_sub(Rd, Rb, BD);
        rvec_sub(Re, Rb, BE);
        rvec_sub(Rd, Rc, CD);
        rvec_sub(Re, Rc, CE);
        rvec_sub(Re, Rd, DE);

        normAB=norm(AB);
        normAC=norm(AC);
        normAD=norm(AD);
        normAE=norm(AE);
        normBC=norm(BC);
        normBD=norm(BD);
        normBE=norm(BE);
        normCD=norm(CD);
        normCE=norm(CE);
        normDE=norm(DE);

        for(i = 0; i < 3; i++)
        {
            if(normAB > epsLS)
                AB[i]/=normAB;
            if(normAC > epsLS)
                AC[i]/=normAC;
            if(normAD > epsLS)
                AD[i]/=normAD;
            if(normAE > epsLS)
                AE[i]/=normAE;
            if(normBC > epsLS)
                BC[i]/=normBC;
            if(normBD > epsLS)
                BD[i]/=normBD;
            if(normBE > epsLS)
                BE[i]/=normBE;
            if(normCD > epsLS)
                CD[i]/=normCD;
            if(normCE > epsLS)
                CE[i]/=normCE;
            if(normDE > epsLS)
                DE[i]/=normDE;
        }

        for ( i = 0; i < nCol*nRow; i++ )
        {
            M[i] = 0.0;
        }

        //Force on particle 1:
        M[nRow*0+0] = AB[0]; M[nRow*1+0] = AC[0]; M[nRow*2+0] = AD[0]; M[nRow*3+0] = AE[0];
        M[nRow*0+1] = AB[1]; M[nRow*1+1] = AC[1]; M[nRow*2+1] = AD[1]; M[nRow*3+1] = AE[1];
        M[nRow*0+2] = AB[2]; M[nRow*1+2] = AC[2]; M[nRow*2+2] = AD[2]; M[nRow*3+2] = AE[2];
        b[0] = Fa[0]; b[1] = Fa[1]; b[2] = Fa[2];

        //Force on particle 2:
        M[nRow*0+3] = -AB[0]; M[nRow*4+3] = BC[0]; M[nRow*5+3] = BD[0]; M[nRow*6+3] = BE[0];
        M[nRow*0+4] = -AB[1]; M[nRow*4+4] = BC[1]; M[nRow*5+4] = BD[1]; M[nRow*6+4] = BE[1];
        M[nRow*0+5] = -AB[2]; M[nRow*4+5] = BC[2]; M[nRow*5+5] = BD[2]; M[nRow*6+5] = BE[2];
        b[3] = Fb[0]; b[4] = Fb[1]; b[5] = Fb[2];

        //Force on particle 3:
        M[nRow*1+6] = -AC[0]; M[nRow*4+6] = -BC[0]; M[nRow*7+6] = CD[0]; M[nRow*8+6] = CE[0];
        M[nRow*1+7] = -AC[1]; M[nRow*4+7] = -BC[1]; M[nRow*7+7] = CD[1]; M[nRow*8+7] = CE[1];
        M[nRow*1+8] = -AC[2]; M[nRow*4+8] = -BC[2]; M[nRow*7+8] = CD[2]; M[nRow*8+8] = CE[2];
        b[6] = Fc[0]; b[7] = Fc[1]; b[8] = Fc[2];

        //Force on particle 4:
        M[nRow*2+9]  = -AD[0]; M[nRow*5+9]  = -BD[0]; M[nRow*7+9]  = -CD[0]; M[nRow*9+9]  = DE[0];
        M[nRow*2+10] = -AD[1]; M[nRow*5+10] = -BD[1]; M[nRow*7+10] = -CD[1]; M[nRow*9+10] = DE[1];
        M[nRow*2+11] = -AD[2]; M[nRow*5+11] = -BD[2]; M[nRow*7+11] = -CD[2]; M[nRow*9+11] = DE[2];
        b[9] = Fd[0]; b[10] = Fd[1]; b[11] = Fd[2];

        //Force on particle 5:
        M[nRow*3+12] = -AE[0]; M[nRow*6+12] = -BE[0]; M[nRow*8+12] = -CE[0]; M[nRow*9+12] = -DE[0];
        M[nRow*3+13] = -AE[1]; M[nRow*6+13] = -BE[1]; M[nRow*8+13] = -CE[1]; M[nRow*9+13] = -DE[1];
        M[nRow*3+14] = -AE[2]; M[nRow*6+14] = -BE[2]; M[nRow*8+14] = -CE[2]; M[nRow*9+14] = -DE[2];
        b[12] = Fe[0]; b[13] = Fe[1]; b[14] = Fe[2];

        /* Query and allocate the optimal workspace */
        lwork = -1;
        F77_FUNC(dgelsd,DGELSD)(&nRow, &nCol, &nRHS, M, &nRow, b, &nRow, s, &eps1, &rank, &wkopt, &lwork, iwork, &info );
        lwork = (int)wkopt;
        work = (real*) malloc( lwork*sizeof(real) );
        /* Solve the equations A*X = B */
        // Least-Squares solution to the system 
        F77_FUNC(dgelsd,DGELSD)(&nRow, &nCol, &nRHS, M, &nRow, b, &nRow, s, &eps1, &rank, work, &lwork, iwork, &info );

        //If cCFD project the least squares CFD to the shape space
        if(grid->fdecomp == encCFD)
        {
            //Calculate the normal to the Shape Space
            gmxLS_ShapeSpace5Normal(normAB,normAC,normAD,normAE,normBC,normBD,normBE,normCD,normCE,normDE,CaleyMengerNormal);

            //Covariant derivative:
            prod = 0.0;
            for ( i = 0; i < nCol; i++ )
            {
                prod += b[i]*CaleyMengerNormal[i];
            }

            for ( i = 0; i < nCol; i++ )
            {
                b[i] = b[i] - prod * CaleyMengerNormal[i];
            }
        }

        // Sum the 10 contributions to the stress
        lab = b[0];
        lac = b[1];
        lad = b[2];
        lae = b[3];
        lbc = b[4];
        lbd = b[5];
        lbe = b[6];
        lcd = b[7];
        lce = b[8];
        lde = b[9];

        Fij[0] = lab * AB[0]; Fij[1] = lab * AB[1]; Fij[2] = lab * AB[2];
        gmxLS_distribute_interaction(grid, Ra, Rb, Fij,0);
        Fij[0] = lac * AC[0]; Fij[1] = lac * AC[1]; Fij[2] = lac * AC[2];
        gmxLS_distribute_interaction(grid, Ra, Rc, Fij,0);
        Fij[0] = lad * AD[0]; Fij[1] = lad * AD[1]; Fij[2] = lad * AD[2];
        gmxLS_distribute_interaction(grid, Ra, Rd, Fij,0);
        Fij[0] = lae * AE[0]; Fij[1] = lae * AE[1]; Fij[2] = lae * AE[2];
        gmxLS_distribute_interaction(grid, Ra, Re, Fij,0);
        Fij[0] = lbc * BC[0]; Fij[1] = lbc * BC[1]; Fij[2] = lbc * BC[2];
        gmxLS_distribute_interaction(grid, Rb, Rc, Fij,0);
        Fij[0] = lbd * BD[0]; Fij[1] = lbd * BD[1]; Fij[2] = lbd * BD[2];
        gmxLS_distribute_interaction(grid, Rb, Rd, Fij,0);
        Fij[0] = lbe * BE[0]; Fij[1] = lbe * BE[1]; Fij[2] = lbe * BE[2];
        gmxLS_distribute_interaction(grid, Rb, Re, Fij,0);
        Fij[0] = lcd * CD[0]; Fij[1] = lcd * CD[1]; Fij[2] = lcd * CD[2];
        gmxLS_distribute_interaction(grid, Rc, Rd, Fij,0);
        Fij[0] = lce * CE[0]; Fij[1] = lce * CE[1]; Fij[2] = lce * CE[2];
        gmxLS_distribute_interaction(grid, Rc, Re, Fij,0);
        Fij[0] = lde * DE[0]; Fij[1] = lde * DE[1]; Fij[2] = lde * DE[2];
        gmxLS_distribute_interaction(grid, Rd, Re, Fij,0);

        free(work);
    }
    else if(grid->fdecomp == enGLD)
    {
        Fij[0] = (Fa[XX]-Fb[XX])/5.0; Fij[1] = (Fa[YY]-Fb[YY])/5.0; Fij[2] = (Fa[ZZ]-Fb[ZZ])/5.0;
        gmxLS_distribute_interaction(grid, Ra, Rb, Fij,0);
        Fij[0] = (Fa[XX]-Fc[XX])/5.0; Fij[1] = (Fa[YY]-Fc[YY])/5.0; Fij[2] = (Fa[ZZ]-Fc[ZZ])/5.0;
        gmxLS_distribute_interaction(grid, Ra, Rc, Fij,0);
        Fij[0] = (Fa[XX]-Fd[XX])/5.0; Fij[1] = (Fa[YY]-Fd[YY])/5.0; Fij[2] = (Fa[ZZ]-Fd[ZZ])/5.0;
        gmxLS_distribute_interaction(grid, Ra, Rd, Fij,0);
        Fij[0] = (Fa[XX]-Fe[XX])/5.0; Fij[1] = (Fa[YY]-Fe[YY])/5.0; Fij[2] = (Fa[ZZ]-Fe[ZZ])/5.0;
        gmxLS_distribute_interaction(grid, Ra, Re, Fij,0);
        Fij[0] = (Fb[XX]-Fc[XX])/5.0; Fij[1] = (Fb[YY]-Fc[YY])/5.0; Fij[2] = (Fb[ZZ]-Fc[ZZ])/5.0;
        gmxLS_distribute_interaction(grid, Rb, Rc, Fij,0);
        Fij[0] = (Fb[XX]-Fd[XX])/5.0; Fij[1] = (Fb[YY]-Fd[YY])/5.0; Fij[2] = (Fb[ZZ]-Fd[ZZ])/5.0;
        gmxLS_distribute_interaction(grid, Rb, Rd, Fij,0);
        Fij[0] = (Fb[XX]-Fe[XX])/5.0; Fij[1] = (Fb[YY]-Fe[YY])/5.0; Fij[2] = (Fb[ZZ]-Fe[ZZ])/5.0;
        gmxLS_distribute_interaction(grid, Rb, Re, Fij,0);
        Fij[0] = (Fc[XX]-Fd[XX])/5.0; Fij[1] = (Fc[YY]-Fd[YY])/5.0; Fij[2] = (Fc[ZZ]-Fd[ZZ])/5.0;
        gmxLS_distribute_interaction(grid, Rc, Rd, Fij,0);
        Fij[0] = (Fc[XX]-Fe[XX])/5.0; Fij[1] = (Fc[YY]-Fe[YY])/5.0; Fij[2] = (Fc[ZZ]-Fe[ZZ])/5.0;
        gmxLS_distribute_interaction(grid, Rc, Re, Fij,0);
        Fij[0] = (Fd[XX]-Fe[XX])/5.0; Fij[1] = (Fd[YY]-Fe[YY])/5.0; Fij[2] = (Fd[ZZ]-Fe[ZZ])/5.0;
        gmxLS_distribute_interaction(grid, Rd, Re, Fij,0);
    }
    else if (grid->fdecomp == enMOP)
    {
        printf("ERROR: 5-particle potentials (CMAP) for MOP haven't been programmed\n");
        exit(1);
    }
}
//----------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------
// AUXILIARY FUNCTIONS

// Modulo operation
int gmxLS_modulo (int a, int b)
{
    int ret;
    ret = a % b;
    if(ret < 0)
    ret+=b;
    return ret;
}

// Finds the indices on the grid for a given set of coordinates
void gmxLS_grid_givecoord(gmxLS_locals_grid_t *grid,rvec pt, int *i, int *j, int *k, int pbc)
{

    int nx,ny,nz;
    real rxx,ryx,ryy,rzx,rzy,rzz;

    /*Looks up indices for grid:
    fractional indices are invbox * coordinate;
    grid indices are then nx*f_ind[XX], etc.
    */
    rxx = grid->invbox[XX][XX];
    ryx = grid->invbox[YY][XX];
    ryy = grid->invbox[YY][YY];
    rzx = grid->invbox[ZZ][XX];
    rzy = grid->invbox[ZZ][YY];
    rzz = grid->invbox[ZZ][ZZ];

    nx = grid->nx;
    ny = grid->ny;
    nz = grid->nz;

    *i = nx * pt[XX] * rxx;
    *j = ny * pt[YY] * ryy;
    *k = nz * pt[ZZ] * rzz;

    if(pbc) {
    *i = gmxLS_modulo(*i,nx);
    *j = gmxLS_modulo(*j,ny);
    *k = gmxLS_modulo(*k,nz);
    }

    if(pt[0] < 0) *i =*i-1;
    if(pt[1] < 0) *j =*j-1;
    if(pt[2] < 0) *k =*k-1;
}

// Duplicate of grid_givecoord that uses the box directly, needed by g_density3D
void gmxLS_grid_givecoord2(matrix box, int nx, int ny, int nz,rvec pt, int *i, int *j, int *k, int pbc)
{

    real rxx,ryx,ryy,rzx,rzy,rzz;

    /*Looks up indices for grid:
    fractional indices are invbox * coordinate;
    grid indices are then nx*f_ind[XX], etc.
    */
    real tmp=1.0/(box[XX][XX]*box[YY][YY]*box[ZZ][ZZ]);
    rxx = box[YY][YY]*box[ZZ][ZZ]*tmp;
    ryy = box[XX][XX]*box[ZZ][ZZ]*tmp;
    rzz = box[XX][XX]*box[YY][YY]*tmp;

    *i = nx * pt[XX] * rxx;
    *j = ny * pt[YY] * ryy;
    *k = nz * pt[ZZ] * rzz;

    if(pbc) {
    *i = gmxLS_modulo(*i,nx);
    *j = gmxLS_modulo(*j,ny);
    *k = gmxLS_modulo(*k,nz);
    }

    if(pt[0] < 0) *i =*i-1;
    if(pt[1] < 0) *j =*j-1;
    if(pt[2] < 0) *k =*k-1;
}
//----------------------------------------------------------------------------------------
// FIVE BODY POTENTIALS -> CALEY-MENGER DERIVATIVES FOR cCFD
// Calculate the derivative of the Caley-Menger determinant for the 5-particles case with respect to d12
real gmxLS_CaleyMenger5Der(real d12,real d13,real d14,real d15,real d23,real d24,real d25,real d34,real d35,real d45)
{
    return -4.0* d12 *( d45*d45*(d35*d35*(-(2.0*d12*d12) + d23*d23 + d24*d24 - 2.0*d34*d34) + d34*d34*(-(2.0*d12*d12) + d23*d23 + d25*d25) + d13*d13*(-(2.0*d23*d23) + d24*d24 + d25*d25 + d34*d34 + d35*d35)) + d45*d45*d45*d45*(-(-d12*d12 + d13*d13 + d23*d23)) - (d34 - d35)*(d34 + d35)*(d34*d34*(d25*d25 - d12*d12) + d35*d35*(d12 - d24)*(d12 + d24) + d13*d13*(d24 - d25)*(d24 + d25)) + d14*d14*(d23*d23*(-d34*d34 + d35*d35 + d45*d45) + d35*d35*(-(2.0*d24*d24) + d34*d34 - d35*d35 + d45*d45) + d25*d25*(d34*d34 + d35*d35 - d45*d45)) + d15*d15*(d23*d23*(d34*d34 - d35*d35 + d45*d45) + d24*d24*(d34*d34 + d35*d35 - d45*d45) + d34*d34*(-(2.0*d25*d25) - d34*d34 + d35*d35 + d45*d45)));
}

//Calculate the normal to the shape space for the 5-particles case
void gmxLS_ShapeSpace5Normal(real d12,real d13,real d14,real d15,real d23,real d24,real d25,real d34,real d35,real d45, real *normal)
{
    real norm;
    int nDist = 10;
    int i;

    real d12_, d13_, d14_, d15_, d23_, d24_, d25_, d34_, d35_, d45_;

    //No change
    d12_ = d12; d13_=d13; d14_=d14; d15_=d15; d23_=d23; d24_=d24; d25_=d25; d34_=d34; d35_=d35; d45_=d45;
    normal[0] = gmxLS_CaleyMenger5Der(d12_,d13_,d14_,d15_,d23_,d24_,d25_,d34_,d35_,d45_);

    d12_ = d13; d13_=d12; d14_=d14; d15_=d15; d23_=d23; d24_=d34; d25_=d35; d34_=d24; d35_=d25; d45_=d45;
    normal[1] = gmxLS_CaleyMenger5Der(d12_,d13_,d14_,d15_,d23_,d24_,d25_,d34_,d35_,d45_);

    d12_ = d14; d13_=d12; d14_=d13; d15_=d15; d23_=d24; d24_=d34; d25_=d45; d34_=d23; d35_=d25; d45_=d35;
    normal[2] = gmxLS_CaleyMenger5Der(d12_,d13_,d14_,d15_,d23_,d24_,d25_,d34_,d35_,d45_);

    d12_ = d15; d13_=d12; d14_=d13; d15_=d14; d23_=d25; d24_=d35; d25_=d45; d34_=d23; d35_=d24; d45_=d34;
    normal[3] = gmxLS_CaleyMenger5Der(d12_,d13_,d14_,d15_,d23_,d24_,d25_,d34_,d35_,d45_);

    d12_ = d23; d13_=d12; d14_=d24; d15_=d25; d23_=d13; d24_=d34; d25_=d35; d34_=d14; d35_=d15; d45_=d45;
    normal[4] = gmxLS_CaleyMenger5Der(d12_,d13_,d14_,d15_,d23_,d24_,d25_,d34_,d35_,d45_);

    d12_ = d24; d13_=d12; d14_=d23; d15_=d25; d23_=d14; d24_=d34; d25_=d45; d34_=d13; d35_=d15; d45_=d35;
    normal[5] = gmxLS_CaleyMenger5Der(d12_,d13_,d14_,d15_,d23_,d24_,d25_,d34_,d35_,d45_);

    d12_ = d25; d13_=d12; d14_=d23; d15_=d24; d23_=d15; d24_=d35; d25_=d45; d34_=d13; d35_=d14; d45_=d34;
    normal[6] = gmxLS_CaleyMenger5Der(d12_,d13_,d14_,d15_,d23_,d24_,d25_,d34_,d35_,d45_);

    d12_ = d34; d13_=d13; d14_=d23; d15_=d35; d23_=d14; d24_=d24; d25_=d45; d34_=d12; d35_=d15; d45_=d25;
    normal[7] = gmxLS_CaleyMenger5Der(d12_,d13_,d14_,d15_,d23_,d24_,d25_,d34_,d35_,d45_);

    d12_ = d35; d13_=d13; d14_=d23; d15_=d34; d23_=d15; d24_=d25; d25_=d45; d34_=d12; d35_=d14; d45_=d24;
    normal[8] = gmxLS_CaleyMenger5Der(d12_,d13_,d14_,d15_,d23_,d24_,d25_,d34_,d35_,d45_);

    d12_ = d45; d13_=d14; d14_=d24; d15_=d34; d23_=d15; d24_=d25; d25_=d35; d34_=d12; d35_=d13; d45_=d23;
    normal[9] = gmxLS_CaleyMenger5Der(d12_,d13_,d14_,d15_,d23_,d24_,d25_,d34_,d35_,d45_);

    norm = 0.0;
    for ( i = 0; i < nDist; i++ )
    {
        norm += normal[i]*normal[i];
    }
    norm = sqrt(norm);

    if (norm > epsLS)
    {
        for ( i = 0; i < nDist; i++ )
        {
            normal[i] /= norm;
        }
    }
}
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
