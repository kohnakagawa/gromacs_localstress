/* 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <math.h>
#include <ctype.h>

#include "sysstuff.h"
#include "string.h"
#include "string2.h"
#include "typedefs.h"
#include "smalloc.h"
#include "macros.h"
#include "gstat.h"
#include "vec.h"
#include "xvgr.h"
#include "pbc.h"
#include "copyrite.h"
#include "futil.h"
#include "statutil.h"
#include "index.h"
#include "tpxio.h"
#include "physics.h"
#include "gmx_ana.h"
#include "localstress.h"

void calc_density3D(const char *fn, atom_id **index, int gnx[], double **slDensity,
                  t_topology *top, int ePBC,int nr_grps, int *gridx, int *gridy, int *gridz,
                  real *gridsp, matrix avgbox, const output_env_t oenv)
{
  rvec *x0;              /* coordinates without pbc */
  matrix box;            /* box (3x3) */
  double invvol,mass,invgridsp,dummy1,dummy2,factor;
  int natoms;            /* nr. atoms in trj */
  int tx,ty,tz,sx,sy,sz,nx,ny,nz,d,dd;
  t_trxstatus *status;  
  int  **slCount,         /* nr. of atoms in one slice for a group */
      i,ii,iii,j,jj,jjj,k,kk,kkk,l,n,               /* loop indices */
      teller = 0,
      nr_frames = 0;     /* number of frames */
  real t,x,y,z;
  char *buf;             /* for tmp. keeping atomname */
  gmx_rmpbc_t  gpbc=NULL;
  rvec gridspv,pt;

  if ((natoms = read_first_x(oenv,&status,fn,&t,&x0,box)) == 0)
    gmx_fatal(FARGS,"Could not read coordinates from statusfile\n");

  gpbc = gmx_rmpbc_init(&top->idef,ePBC,top->atoms.nr,box);
  for (d = 0; d < 3; d++)
      for (dd = 0; dd < 3; dd++)
        avgbox[d][dd] = 0.0;
  if (*gridx == 0)
    *gridx = box[XX][XX]/(*gridsp);
  if (*gridy == 0)
    *gridy = box[YY][YY]/(*gridsp);
  if (*gridz == 0)
    *gridz = box[ZZ][ZZ]/(*gridsp);
    
  printf("gridx = %d; gridy = %d, gridz = %d \n", *gridx, *gridy, *gridz);
  nx = *gridx;
  ny = *gridy;
  nz = *gridz;

  snew(*slDensity, nx*ny*nz);
  /*********** Start processing trajectory ***********/
  do {
    gmx_rmpbc(gpbc,natoms,box,x0);

    for (d = 0; d < 3; d++)
      for (dd = 0; dd < 3; dd++)
        avgbox[d][dd] += box[d][dd];

    invvol = nx*ny*nz/(box[XX][XX]*box[YY][YY]*box[ZZ][ZZ]);
    teller++;
    
    // Define the grid spacing in each direction
    gridspv[XX] = box[XX][XX]/nx;
    gridspv[YY] = box[YY][YY]/ny;
    gridspv[ZZ] = box[ZZ][ZZ]/nz;
    invgridsp = 1/(gridspv[XX]*gridspv[YY]*gridspv[ZZ]);

    for (n = 0; n < nr_grps; n++) {
      for (l = 0; l < gnx[n]; l++) {   /* loop over all atoms in index file */
        pt[0] = x0[index[n][l]][0];
        pt[1] = x0[index[n][l]][1];
        pt[2] = x0[index[n][l]][2];
        
        // Get the coordinates of the point in the grid
        gmxLS_grid_givecoord2(box,nx,ny,nz,pt,&ii,&jj,&kk,0);
        iii=ii;
        jjj=jj;
        kkk=kk;
        
        // Spread it
        for(i=1;i>=-1;i-=2){
          iii+=i;
          dummy1 = i * invgridsp * (pt[0]-(ii+0.5*(1-i))*gridspv[0]);
          for(j=1;j>=-1;j-=2){
            jjj+=j;
            dummy2 = dummy1 * j * (pt[1]-(jj+0.5*(1-j))*gridspv[1]);
            for(k=1;k>=-1;k-=2){
              kkk+=k;
              factor = dummy2 * k * (pt[2]-(kk+0.5*(1-k))*gridspv[2]);
              mass = top->atoms.atom[index[n][l]].m;
              (*slDensity)[gmxLS_modulo(iii,nx)*nz*ny+gmxLS_modulo(jjj,ny)*nz+gmxLS_modulo(kkk,nz)] += factor*mass*invvol;
            }
          }
        }
      }
    }
    nr_frames++;
  } while (read_next_x(oenv,status,&t,natoms,x0,box));
  gmx_rmpbc_done(gpbc);

  for (d = 0; d < 3; d++)
    for (dd = 0; dd < 3; dd++)
      avgbox[d][dd] /= nr_frames;

  /*********** done with status file **********/
  close_trj(status);

  fprintf(stderr,"\nRead %d frames from trajectory. Calculating density\n", nr_frames);
  for (i = 0; i < nx*ny*nz; i++)
    (*slDensity)[i] /= nr_frames;

  sfree(x0);  /* free memory used by coordinate array */
}

int gmx_density3D(int argc,char *argv[])
{
  const char *desc[] = {
    "Computes the density in 3D over the simulation volume and stores it"
    "in a grid with a given spacing. Use the bindens2ncdf program in the"
    "src/tensortools folder to convert the output into a netcdf file"
    "that can be read by a number of programs such as paraview."
  };

  output_env_t oenv;
  static const char *dens_opt[] = 
    { NULL, "mass", "number", "charge", NULL };
  static int  ngrps   = 1;       /* nr. of groups              */
  real gridsp;           /* grid spacing in nm         */
  int gridx,gridy,gridz;  /* grid size in x,y, and z    */
  gridx = 0;
  gridy = 0;
  gridz = 0;
  gridsp = 0.05;
  t_pargs pa[] = {
    { "-gridsp",  FALSE, etREAL, {&gridsp},
      "Divide the box in 3D using a grid with a spacing given by gridsp in nm" },
    { "-gridx",  FALSE, etINT, {&gridx},
      "Override the size of the grid in the x direction" },
    { "-gridy",  FALSE, etINT, {&gridy},
      "Override the size of the grid in the y direction" },
    { "-gridz",  FALSE, etINT, {&gridz},
      "Override the size of the grid in the z direction" },
    { "-dens",    FALSE, etENUM, {dens_opt},
      "Density"}
  };

  const char *bugs[] = {
    "When calculating electron densities, atomnames are used instead of types. This is bad.",
  };
  
  double *density;      /* density array         */
  
  char **grpname;        /* groupnames                 */
  int  *ngx;             /* sizes of groups            */
  t_topology *top;       /* topology                   */
  int  ePBC, d, dd;
  atom_id   **index;     /* indices for all groups     */
  int  i,bdouble;
  bdouble = 1;
  double dbox[3][3];
  matrix avgbox;
  FILE *outfp;

  t_filenm  fnm[] = {    /* files for g_density        */
    { efTRX, "-f", NULL,  ffREAD },
    { efNDX, NULL, NULL,  ffOPTRD },
    { efTPX, NULL, NULL,  ffREAD },
    { efDAT,"-o","density",ffWRITE }
  };

#define NFILE asize(fnm)
  
  CopyRight(stderr,argv[0]);

  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME | PCA_BE_NICE,
                    NFILE,fnm,asize(pa),pa,asize(desc),desc,asize(bugs),bugs,
                    &oenv);
  top = read_top(ftp2fn(efTPX,NFILE,fnm),&ePBC);     /* read topology file */

  
  if (dens_opt[0][0] == 'm') {
    for(i=0; (i<top->atoms.nr); i++)
      top->atoms.atom[i].m = top->atoms.atom[i].m*AMU/(NANO*NANO*NANO);
  } else if (dens_opt[0][0] == 'n') {
    for(i=0; (i<top->atoms.nr); i++)
      top->atoms.atom[i].m = 1;
  } else if (dens_opt[0][0] == 'c') {
    for(i=0; (i<top->atoms.nr); i++)
      top->atoms.atom[i].m = top->atoms.atom[i].q;
  }

  snew(grpname,ngrps);
  snew(index,ngrps);
  snew(ngx,ngrps);

  get_index(&top->atoms,ftp2fn_null(efNDX,NFILE,fnm),ngrps,ngx,index,grpname);
  calc_density3D(ftp2fn(efTRX,NFILE,fnm),index, ngx, &density, top, ePBC, ngrps, &gridx, &gridy, &gridz, &gridsp, avgbox, oenv);
  
  for (d = 0; d < 3; d++)
      for (dd = 0; dd < 3; dd++)
        dbox[d][dd] = avgbox[d][dd];
  outfp = fopen(ftp2fn(efDAT,NFILE,fnm),"w");
  fwrite(&bdouble,sizeof(int),1,outfp);
  fwrite(dbox,sizeof(double),9,outfp);
  fwrite(&gridx,sizeof(int),1,outfp);
  fwrite(&gridy,sizeof(int),1,outfp);
  fwrite(&gridz,sizeof(int),1,outfp);
  fwrite(density,sizeof(double),gridx*gridy*gridz,outfp);
  fclose(outfp);
  do_view(oenv,opt2fn("-o",NFILE,fnm), "-nxy");       /* view xvgr file */
  thanx(stderr);
  return 0;
}

/*
int main(int argc, char *argv[])
{
  gmx_density(argc,argv);
  return 0;
}
*/
