/*  -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 * 
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
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "typedefs.h"
#include "macros.h"
#include "copyrite.h"
#include "main.h"
#include "statutil.h"
#include "smalloc.h"
#include "futil.h"
#include "smalloc.h"
#include "edsam.h"
#include "mdrun.h"
#include "xmdrun.h"
#include "checkpoint.h"
#include <string.h>
#ifdef GMX_THREADS
#include "thread_mpi.h"
#endif

/* afm stuf */
#include "pull.h"

int main(int argc,char *argv[])
{

  const char *desc[] = {
    "\n"
    "--------------------------------------------------------------------------------\n"
    "Custom GROMACS version to compute the local stress tensor in 3D over a\n"
    "simulation. This program computes the local stress based on the Hardy\n"
    "stress definition (see Ch.8 in Tadmor & Miller (2011) Modeling Materials:\n"
    "Continuum, Atomistic and Multiscale Techniques). It is also possible to \n"
    "compute the virial stress per atom (J. Chem. Phys. 131, 154107 (2009)). We \n"
    "have patched the GROMACS source code (v4.5.5) at various locations in \n"
    "the routines that calculate particles forces and velocities. Our patches rely\n"
    "heavily on the code structure and functions implemented in a previous local\n"
    "stress code (obtained from\n"
    "http://repo.or.cz/w/gromacs.git/shortlog/refs/heads/release-4-5-localstress\n"
    "and ftp://ftp.gromacs.org/pub/tmp/gromacs-4.0.2_localstress.tar.gz). This\n"
    "previous code implements the methodology outlined in\n"
    "\n"
    "S. Ollila et al. Phys. Rev. Lett. 102, 078101 (2009). 3D Pressure Field in\n"
    "Lipid Membranes and Membrane-Protein Complexes.\n"
    "\n"
    "In the current custom GROMACS version we implement the methodology outlined\n"
    "in\n"
    "\n"
    "J. M. Vanegas, A. Torres-Sanchez, and M. Arroyo, J. Chem. Theor.\n"
    "Comput. 10, 691-702 (2014). Importance of Force Decomposition for Local\n"
    "Stress Calculations in Biomembrane Molecular Simulations.\n"
    "and in\n"
    "A. Torres-Sanchez, J. M. Vanegas, M. Arroyo, submitted to PRL (2015)\n"
    "\n"
    "The main differences of our implementation compared the previous one\n"
    "include:\n"
    "\n"
    "- Decomposition of multi-body potential forces using:\n"
    "      1. Covariant central force decomposition (cCFD),\n"
    "         Torres-Sanchez A. Vanegas, J. and Arroyo, M. (Submitted),\n" 
    "      2. Non-covariant central force decomposition (nCFD), \n"
    "         N. C. Admal and E. B. Tadmor; J. Elast. 100, 63-143, 2010\n"
    "      3. Goetz-Lipowsky decomposition (GLD)\n"
    "         R. Goetz and R. J. Lipowsky; J. Chem. Phys. 108, 7397-7409, 1998.\n"
    "      4. Decomposition on geometric centers (GMC)\n"
    "         H. Heinz; W. Paul; K. Binder; Phys. Rev. E. 72 066704 (2005)\n"
    "The choice of decomposition can produce drastically different stress tensors\n"
    "due to the contributions from multibody potentials. Any flavour of the CFD\n"
    "always results in a symmetric stress tensor by definition (consistent with\n"
    "the continuum conservation of angular momentum), while GLD in general does\n"
    "not. See the text by Vanegas et al. and Torres-Sanchez et al. for more\n"
    "details. cCFD and nCFD differ in the treatment of 5-body body potentials\n"
    "-CMAP- (see Torres-Sanchez et al.). The decomposition on geometry centers\n"
    "does not satisfy balance of linear momentum and angular momentum (see Torres\n"
    " -Sanchez et al.) and therefore its use is highly discouraged.\n"
    "WE "
    "\n"
    "- Ability to output the total or individual contributions to the local\n"
    "stress such as those from vdw, electrostatics, angles, and others.\n"
    "\n"
    "- Finite discretization over a rectangular grid using tri-linear weight\n"
    "functions, which result in smoother stress fields and also make the\n"
    "discretization exact regardless of the grid size.\n"
    "\n"
    "- Consistent treatment of forces arising from bond constraints (LINCS,\n"
    "SETTLE and SHAKE algorithms).\n"
    "\n"
    "- Virial stress per atom. For comparison, we have recently included the \n"
    "the virial stress per atom as defined in \n\n"
    "A. P. Thompson et al., J. Chem. Phys. 131, 154107 (2009)\n\n"
    "The computation of the virial stress per atom is much faster than the\n" 
    "computation of the IKN stress in a grid. However, we strongly recommend\n"
    "using the virial stress per atom only for rapid visualization, as it\n" 
    "does not satisfy balance of linear momentum (see Torres-Sanchez et al.).\n"
    "\n"
    "If you publish results using this code, we kindly ask you to cite both the\n"
    "paper by Ollila et al. and by Vanegas et al.\n"
    "\n"
    "Brief usage (see Local_stress.pdf file for more details):\n"
    "\n"
    "This program is not meant to run simulations, but to reanalyze a trajectory\n"
    "using the -rerun option. You need to have saved both the positions and\n"
    "velocities together in order for the analysis to work correctly. The output\n"
    "files created by this program (.dat0) are binary files which can be\n"
    "processed with the utilities in the tensor_tools folder in the program\n"
    "source.\n"
    "\n"
    "Remember that this program is very slow compared to regular gromacs as it\n"
    "is not optimized. Also, depending on the size of your system, it can use a\n"
    "large amount of memory. The memory required depends on your grid size, so\n"
    "if youhave a 10x10x10nm system and you choose a grid spacing of 0.1 nm (the\n"
    "default), then you have a grid with 1,000,000 elements. Each one of those\n"
    "elements needs 72 bytes to store the 3x3 tensor in double precision, and we\n"
    "have 2 grids in a run, which adds up to about 140 MB.\n"
    "\n"
    "Although this code cannot be run in parallel directly, each trajectory\n"
    "frame can be analyzed completely independently of every other frame, so you\n"
    "can split your trajectory into separate files (with the same number of\n"
    "frames) and analyze each file separately. You can then average the output\n"
    ".dat0 files using the aver3Dbin utility.\n"
    "\n"
    "The code can output the individual (e.g. vdw, coulomb, angles, etc.) or\n"
    "total contributions to the stress tensor with the -lscont flag. The grid\n"
    "size in every box dimension can be independently controlled using the\n"
    "-lsgridx, -lsgridy, and -lsgridz flags, or by simply specifying the grid\n"
    "spacing with the -localsgrid flag. You can also combine these flags\n"
    "together to optimize the grid, so for example if you are only interested in\n"
    "the stress profile along a given dimension (e.g. z), then you can use the\n"
    "combination of -lsgridx 1 -lsgridy 1 -localsgrid 0.1 to obtain a 1D grid\n"
    "with a spacing of 0.1 nm along the z axis. The choice of force\n"
    "decomposition can be selected with the -lsfd flag. We recommend using the\n"
    "covariant central force decomposition (cCFD). You can also activate the \n"
    "computation of the stress per atom, for which no force decomposition is\n"
    "required.\n"
    "\n"
    "If you have any questions, comments, or bugs please send an email to:\n"
    "juan.m.vanegas@gmail.com or torres.sanchez.a@gmail.com\n"
    "Juan M. Vanegas and Alejandro Torres-Sanchez - 2015\n"
    "--------------------------------------------------------------------------------\n"
  };
  t_commrec    *cr;
  t_filenm fnm[] = {
    { efTPX, NULL,      NULL,       ffREAD },
    { efTRN, "-o",      NULL,       ffWRITE },
    { efXTC, "-x",      NULL,       ffOPTWR },
    { efCPT, "-cpi",    NULL,       ffOPTRD },
    { efCPT, "-cpo",    NULL,       ffOPTWR },
    { efSTO, "-c",      "confout",  ffWRITE },
    { efEDR, "-e",      "ener",     ffWRITE },
    { efLOG, "-g",      "md",       ffWRITE },
    { efXVG, "-dhdl",   "dhdl",     ffOPTWR },
    { efXVG, "-field",  "field",    ffOPTWR },
    { efXVG, "-table",  "table",    ffOPTRD },
    { efXVG, "-tablep", "tablep",   ffOPTRD },
    { efXVG, "-tableb", "table",    ffOPTRD },
    { efTRX, "-rerun",  "rerun",    ffOPTRD },
    { efXVG, "-tpi",    "tpi",      ffOPTWR },
    { efXVG, "-tpid",   "tpidist",  ffOPTWR },
    { efEDI, "-ei",     "sam",      ffOPTRD },
    { efEDO, "-eo",     "sam",      ffOPTWR },
    { efGCT, "-j",      "wham",     ffOPTRD },
    { efGCT, "-jo",     "bam",      ffOPTWR },
    { efXVG, "-ffout",  "gct",      ffOPTWR },
    { efXVG, "-devout", "deviatie", ffOPTWR },
    { efXVG, "-runav",  "runaver",  ffOPTWR },
    { efXVG, "-px",     "pullx",    ffOPTWR },
    { efXVG, "-pf",     "pullf",    ffOPTWR },
    { efMTX, "-mtx",    "nm",       ffOPTWR },
    { efNDX, "-dn",     "dipole",   ffOPTWR },
    { efDAT, "-ols",    "localstress", ffWRITE },
    { efRND, "-multidir",NULL,      ffOPTRDMULT}
  };
#define NFILE asize(fnm)

  /* Command line options ! */
  gmx_bool bCart        = FALSE;
  gmx_bool bPPPME       = FALSE;
  gmx_bool bPartDec     = FALSE;
  gmx_bool bDDBondCheck = TRUE;
  gmx_bool bDDBondComm  = TRUE;
  gmx_bool bVerbose     = FALSE;
  gmx_bool bCompact     = TRUE;
  gmx_bool bSepPot      = FALSE;
  gmx_bool bRerunVSite  = FALSE;
  gmx_bool bIonize      = FALSE;
  gmx_bool bConfout     = TRUE;
  gmx_bool bReproducible = FALSE;
    
  int  npme=-1;
  int  nmultisim=0;
  int  nstglobalcomm=-1;
  int  repl_ex_nst=0;
  int  repl_ex_seed=-1;
  int  nstepout=100;
  int  nthreads=0; /* set to determine # of threads automatically */
  int  resetstep=-1;
  
  rvec realddxyz={0,0,0};
  const char *ddno_opt[ddnoNR+1] =
    { NULL, "interleave", "pp_pme", "cartesian", NULL };
  const char *dddlb_opt[] =
    { NULL, "auto", "no", "yes", NULL };
  real rdd=0.0,rconstr=0.0,dlb_scale=0.8,pforce=-1;
  char *ddcsx=NULL,*ddcsy=NULL,*ddcsz=NULL;
  real cpt_period=15.0,max_hours=-1;
  gmx_bool bAppendFiles=TRUE;
  gmx_bool bKeepAndNumCPT=FALSE;
  gmx_bool bResetCountersHalfWay=FALSE;
  output_env_t oenv=NULL;
  const char *deviceOptions = "";
  real localsgridspacing=0.1;
  int nstlocals=0;
  int localsgridx=0;
  int localsgridy=0;
  int localsgridz=0;
  //enum { enSel, enAll, enVdw, enCoul, enAngles, enBonds, enDihp, enDihi, enDihrb, enLincs, enSettle, enShake, enVel, enNR };
  //const char *localsenum[enNR+1] =
  //    { NULL, "all", "vdw", "coul", "angles", "bonds", "dihp", "dihi", "dihrb", "lincs", "settle", "shake", "vel", NULL};
  char *localsenum   = "all";
  char *localsfdenum = "ccfd";
  char *localssanum  = "spat";
  int compareLimit = 5;

  t_pargs pa[] = {

    { "-pd",      FALSE, etBOOL,{&bPartDec},
      "HIDDENUse particle decompostion" },
    { "-dd",      FALSE, etRVEC,{&realddxyz},
      "HIDDENDomain decomposition grid, 0 is optimize" },
#ifdef GMX_THREADS
    { "-nt",      FALSE, etINT, {&nthreads},
      "HIDDENNumber of threads to start (0 is guess)" },
#endif
    { "-npme",    FALSE, etINT, {&npme},
      "HIDDENNumber of separate nodes to be used for PME, -1 is guess" },
    { "-ddorder", FALSE, etENUM, {ddno_opt},
      "HIDDENDD node order" },
    { "-ddcheck", FALSE, etBOOL, {&bDDBondCheck},
      "HIDDENCheck for all bonded interactions with DD" },
    { "-ddbondcomm", FALSE, etBOOL, {&bDDBondComm},
      "HIDDENUse special bonded atom communication when [TT]-rdd[tt] > cut-off" },
    { "-rdd",     FALSE, etREAL, {&rdd},
      "HIDDENThe maximum distance for bonded interactions with DD (nm), 0 is determine from initial coordinates" },
    { "-rcon",    FALSE, etREAL, {&rconstr},
      "HIDDENMaximum distance for P-LINCS (nm), 0 is estimate" },
    { "-dlb",     FALSE, etENUM, {dddlb_opt},
      "HIDDENDynamic load balancing (with DD)" },
    { "-dds",     FALSE, etREAL, {&dlb_scale},
      "HIDDENMinimum allowed dlb scaling of the DD cell size" },
    { "-ddcsx",   FALSE, etSTR, {&ddcsx},
      "HIDDENThe DD cell sizes in x" },
    { "-ddcsy",   FALSE, etSTR, {&ddcsy},
      "HIDDENThe DD cell sizes in y" },
    { "-ddcsz",   FALSE, etSTR, {&ddcsz},
      "HIDDENThe DD cell sizes in z" },
    { "-gcom",    FALSE, etINT,{&nstglobalcomm},
      "HIDDENGlobal communication frequency" },
    { "-v",       FALSE, etBOOL,{&bVerbose},  
      "Be loud and noisy" },
    { "-compact", FALSE, etBOOL,{&bCompact},  
      "HIDDENWrite a compact log file" },
    { "-seppot",  FALSE, etBOOL, {&bSepPot},
      "HIDDENWrite separate V and dVdl terms for each interaction type and node to the log file(s)" },
    { "-pforce",  FALSE, etREAL, {&pforce},
      "HIDDENPrint all forces larger than this (kJ/mol nm)" },
    { "-reprod",  FALSE, etBOOL,{&bReproducible},  
      "HIDDENTry to avoid optimizations that affect binary reproducibility" },
    { "-cpt",     FALSE, etREAL, {&cpt_period},
      "HIDDENCheckpoint interval (minutes)" },
    { "-cpnum",   FALSE, etBOOL, {&bKeepAndNumCPT},
      "HIDDENKeep and number checkpoint files" },
    { "-append",  FALSE, etBOOL, {&bAppendFiles},
      "HIDDENAppend to previous output files when continuing from checkpoint instead of adding the simulation part number to all file names" },
    { "-maxh",   FALSE, etREAL, {&max_hours},
      "HIDDENTerminate after 0.99 times this time (hours)" },
    { "-multi",   FALSE, etINT,{&nmultisim}, 
      "HIDDENDo multiple simulations in parallel" },
    { "-replex",  FALSE, etINT, {&repl_ex_nst}, 
      "HIDDENAttempt replica exchange every # steps" },
    { "-reseed",  FALSE, etINT, {&repl_ex_seed}, 
      "HIDDENSeed for replica exchange, -1 is generate a seed" },
    { "-localsgrid",  FALSE, etREAL, {&localsgridspacing},
      "Spacing for local stress grid (default = 0.1 nm)" },
    { "-nstlp",  FALSE, etINT, {&nstlocals},
      "HIDDENFrequency of writing local stress grid to file (default = 0)" },
    { "-lsgridx", FALSE, etINT, {&localsgridx},
      "Set the local stress grid size in the x direction (default use box[XX][XX]/localsgrid)"},
    { "-lsgridy", FALSE, etINT, {&localsgridy},
      "Set the local stress grid size in the y direction (default use box[YY][YY]/localsgrid)"},
    { "-lsgridz", FALSE, etINT, {&localsgridz},
      "Set the local stress grid size in the z direction (default use box[ZZ][ZZ]/localsgrid)"},
    { "-lscont", FALSE, etSTR, {&localsenum},
      "Select which contribution to write to output (default = all): all, vdw, coul, angles, bonds, dihp, dihi, dihrb, lincs, settle, shake, cmap, vel"},
    { "-lsfd", FALSE, etSTR, {&localsfdenum},
      "Select the type of force decomposition to be used: ccfd (covariant central force decomposition, default), ncfd (non-covariant central force decomposition), gld (Goetz-Lipowsky decomposition), or gmc (decomposition on geometric centers)"},
    { "-lssa", FALSE, etSTR, {&localssanum},
      "Select the type of stress to calculate: spat (spatial stress from IKN theory, default), atom (stress per atom)"},
    { "-rerunvsite", FALSE, etBOOL, {&bRerunVSite},
      "HIDDENRecalculate virtual site coordinates with [TT]-rerun[tt]" },
    { "-ionize",  FALSE, etBOOL,{&bIonize},
      "HIDDENDo a simulation including the effect of an X-Ray bombardment on your system" },
    { "-confout", FALSE, etBOOL, {&bConfout},
      "HIDDENWrite the last configuration with [TT]-c[tt] and force checkpointing at the last step" },
    { "-stepout", FALSE, etINT, {&nstepout},
      "HIDDENFrequency of writing the remaining runtime" },
    { "-resetstep", FALSE, etINT, {&resetstep},
      "HIDDENReset cycle counters after these many time steps" },
    { "-resethway", FALSE, etBOOL, {&bResetCountersHalfWay},
      "HIDDENReset the cycle counters after half the number of steps or halfway [TT]-maxh[tt]" }
#ifdef GMX_OPENMM
    ,
    { "-device",  FALSE, etSTR, {&deviceOptions},
      "Device option string" }
#endif
  };
  gmx_edsam_t  ed;
  unsigned long Flags, PCA_Flags;
  ivec     ddxyz;
  int      dd_node_order;
  int      localscontrib;
  int      localsfdecomp;
  int      localsspatialatom;
  gmx_bool     bAddPart;
  FILE     *fplog,*fptest;
  int      sim_part,sim_part_fn;
  const char *part_suffix=".part";
  char     suffix[STRLEN];
  int      rc;
  char **multidir=NULL;
  
  

  
  cr = init_par(&argc,&argv);

  if (MASTER(cr))
    CopyRight(stderr, argv[0]);

  PCA_Flags = (PCA_KEEP_ARGS | PCA_NOEXIT_ON_ARGS | PCA_CAN_SET_DEFFNM
	       | (MASTER(cr) ? 0 : PCA_QUIET));
  

  /* Comment this in to do fexist calls only on master
   * works not with rerun or tables at the moment
   * also comment out the version of init_forcerec in md.c 
   * with NULL instead of opt2fn
   */
  /*
     if (!MASTER(cr))
     {
     PCA_Flags |= PCA_NOT_READ_NODE;
     }
     */

  parse_common_args(&argc,argv,PCA_Flags, NFILE,fnm,asize(pa),pa,
                    asize(desc),desc,0,NULL, &oenv);



  /* we set these early because they might be used in init_multisystem() 
     Note that there is the potential for npme>nnodes until the number of
     threads is set later on, if there's thread parallelization. That shouldn't
     lead to problems. */ 
  dd_node_order = nenum(ddno_opt);
  cr->npmenodes = npme;
  

  if (strncmp(localssanum,"spat",4) == 0) {
    localsspatialatom = enSpat;
    printf("\nSelected spatial stress tensor\n");
  }else if (strncmp(localssanum,"atom",4) == 0) {
    localsspatialatom = enAtom;
    printf("\nSelected stress tensor by atom. Will not use force decomposition flag.\n");
  }else{
    printf("\nOption not recognized, will use spatial stress tensor\n");
    localsspatialatom = enSpat;
  }
  
  if (strncmp(localsfdenum,"ccfd",4) == 0) {
    localsfdecomp = encCFD;
    printf("\nSelected force decomposition: %s\n", localsfdenum);
  }else if (strncmp(localsfdenum,"ncfd",4) == 0) {
    localsfdecomp = enCFD;
    printf("\nSelected force decomposition: %s\n", localsfdenum);
  }else if(strncmp(localsfdenum,"gld",4) == 0){
    localsfdecomp = enGLD;
    printf("\nSelected force decomposition: %s\n", localsfdenum);
  }else if(strncmp(localsfdenum,"gmc",4) == 0){
    localsfdecomp = enGMC;
    printf("\nSelected force decomposition: %s\n", localsfdenum);
  }else{
    printf("\nOption not recognized, will use covariant central force decomposition\n");
    localsfdecomp = encCFD;
  }
  printf("\nSelected contribution: %s\n",localsenum);
  if (strncmp(localsenum,"all",5) == 0) {
    printf("\nWill write all contributions to the local stress\n");
    localscontrib = enAll;
  }else if(strncmp(localsenum,"vdw",5) == 0){
    printf("\nWill only write vdw contributions to the local stress\n");
    localscontrib = enVdw;
  }else if(strncmp(localsenum,"coul",5) == 0){
    printf("\nWill only write coulomb contributions to the local stress\n");
    localscontrib = enCoul;
  }else if(strncmp(localsenum,"angles",5) == 0){
    printf("\nWill only write angle contributions to the local stress\n");
    localscontrib = enAngles;
  }else if(strncmp(localsenum,"bonds",5) == 0){
    printf("\nWill only write bonding contributions to the local stress\n");
    localscontrib = enBonds;
  }else if(strncmp(localsenum,"dihp",5) == 0){
    printf("\nWill only write proper dihedral contributions to the local stress\n");
    localscontrib = enDihp;
  }else if(strncmp(localsenum,"dihi",5) == 0){
    printf("\nWill only write inproper dihedral contributions to the local stress\n");
    localscontrib = enDihi;
  }else if(strncmp(localsenum,"dihrb",5) == 0){
    printf("\nWill only write RB dihedral contributions to the local stress\n");
    localscontrib = enDihrb;
  }else if(strncmp(localsenum,"lincs",5) == 0){
    printf("\nWill only write LINCS constraints contributions to the local stress\n");
    localscontrib = enLincs;
  }else if(strncmp(localsenum,"settle",5) == 0){
    printf("\nWill only write SETTLE water constraints contributions to the local stress\n");
    localscontrib = enSettle;
  }else if(strncmp(localsenum,"shake",5) == 0){
    printf("\nWill only write SHAKE constraints contributions to the local stress\n");
    localscontrib = enShake;
  }else if(strncmp(localsenum,"vel",5) == 0){
    printf("\nWill only write velocity contributions to the local stress\n");
    localscontrib = enVel;
  }else if(strncmp(localsenum,"cmap",5) == 0){
    printf("\nWill only write CMAP contributions to the local stress\n");
    localscontrib = enCMAP;
  }else{
    printf("\nOption not recognized, will write all contributions to the local stress\n");
    localscontrib = enAll;
  }

#ifndef GMX_THREADS
  nthreads=1;
#endif

  /* now check the -multi and -multidir option */
  if (opt2bSet("-multidir", NFILE, fnm))
  {
      int i;
      if (nmultisim > 0)
      {
          gmx_fatal(FARGS, "mdrun -multi and -multidir options are mutually exclusive.");
      }
      nmultisim = opt2fns(&multidir, "-multidir", NFILE, fnm);
  }


  if (repl_ex_nst != 0 && nmultisim < 2)
      gmx_fatal(FARGS,"Need at least two replicas for replica exchange (option -multi)");

  if (nmultisim > 1) {
#ifndef GMX_THREADS
    gmx_bool bParFn = (multidir == NULL);
    init_multisystem(cr, nmultisim, multidir, NFILE, fnm, bParFn);
#else
    gmx_fatal(FARGS,"mdrun -multi is not supported with the thread library.Please compile GROMACS with MPI support");
#endif
  }

  bAddPart = !bAppendFiles;

  /* Check if there is ANY checkpoint file available */	
  sim_part    = 1;
  sim_part_fn = sim_part;
  if (opt2bSet("-cpi",NFILE,fnm))
  {
      if (bSepPot && bAppendFiles)
      {
          gmx_fatal(FARGS,"Output file appending is not supported with -seppot");
      }

      bAppendFiles =
                read_checkpoint_simulation_part(opt2fn_master("-cpi", NFILE,
                                                              fnm,cr),
                                                &sim_part_fn,NULL,cr,
                                                bAppendFiles,NFILE,fnm,
                                                part_suffix,&bAddPart);
      if (sim_part_fn==0 && MASTER(cr))
      {
          fprintf(stdout,"No previous checkpoint file present, assuming this is a new run.\n");
      }
      else
      {
          sim_part = sim_part_fn + 1;
      }

      if (MULTISIM(cr))
      {
          check_multi_int(stdout,cr->ms,sim_part,"simulation part");
      }
  } 
  else
  {
      bAppendFiles = FALSE;
  }

  if (!bAppendFiles)
  {
      sim_part_fn = sim_part;
  }

  if (bAddPart)
  {
      /* Rename all output files (except checkpoint files) */
      /* create new part name first (zero-filled) */
      sprintf(suffix,"%s%04d",part_suffix,sim_part_fn);

      add_suffix_to_output_names(fnm,NFILE,suffix);
      if (MASTER(cr))
      {
          fprintf(stdout,"Checkpoint file is from part %d, new output files will be suffixed '%s'.\n",sim_part-1,suffix);
      }
  }

  Flags = opt2bSet("-rerun",NFILE,fnm) ? MD_RERUN : 0;
  Flags = Flags | (bSepPot       ? MD_SEPPOT       : 0);
  Flags = Flags | (bIonize       ? MD_IONIZE       : 0);
  Flags = Flags | (bPartDec      ? MD_PARTDEC      : 0);
  Flags = Flags | (bDDBondCheck  ? MD_DDBONDCHECK  : 0);
  Flags = Flags | (bDDBondComm   ? MD_DDBONDCOMM   : 0);
  Flags = Flags | (bConfout      ? MD_CONFOUT      : 0);
  Flags = Flags | (bRerunVSite   ? MD_RERUN_VSITE  : 0);
  Flags = Flags | (bReproducible ? MD_REPRODUCIBLE : 0);
  Flags = Flags | (bAppendFiles  ? MD_APPENDFILES  : 0); 
  Flags = Flags | (bKeepAndNumCPT ? MD_KEEPANDNUMCPT : 0); 
  Flags = Flags | (sim_part>1    ? MD_STARTFROMCPT : 0); 
  Flags = Flags | (bResetCountersHalfWay ? MD_RESETCOUNTERSHALFWAY : 0);


  /* We postpone opening the log file if we are appending, so we can 
     first truncate the old log file and append to the correct position 
     there instead.  */
  if ((MASTER(cr) || bSepPot) && !bAppendFiles) 
  {
      gmx_log_open(ftp2fn(efLOG,NFILE,fnm),cr,!bSepPot,Flags,&fplog);
      CopyRight(fplog,argv[0]);
      please_cite(fplog,"Hess2008b");
      please_cite(fplog,"Spoel2005a");
      please_cite(fplog,"Lindahl2001a");
      please_cite(fplog,"Berendsen95a");
  }
  else if (!MASTER(cr) && bSepPot)
  {
      gmx_log_open(ftp2fn(efLOG,NFILE,fnm),cr,!bSepPot,Flags,&fplog);
  }
  else
  {
      fplog = NULL;
  }

  ddxyz[XX] = (int)(realddxyz[XX] + 0.5);
  ddxyz[YY] = (int)(realddxyz[YY] + 0.5);
  ddxyz[ZZ] = (int)(realddxyz[ZZ] + 0.5);
  rc = mdrunner(nthreads, fplog,cr,NFILE,fnm,oenv,bVerbose,bCompact,
                nstglobalcomm, ddxyz,dd_node_order,rdd,rconstr,
                dddlb_opt[0],dlb_scale,ddcsx,ddcsy,ddcsz,
                nstepout,resetstep,nmultisim,repl_ex_nst,repl_ex_seed,
                pforce, cpt_period,max_hours,deviceOptions,localsgridspacing,nstlocals,localsgridx,
                localsgridy,localsgridz,localscontrib,localsfdecomp,localsspatialatom,Flags);

  if (gmx_parallel_env_initialized())
      gmx_finalize();

  if (MULTIMASTER(cr)) {
      thanx(stderr);
  }

  /* Log file has to be closed in mdrunner if we are appending to it 
     (fplog not set here) */
  if (MASTER(cr) && !bAppendFiles) 
  {
      gmx_log_close(fplog);
  }

  return rc;
}

