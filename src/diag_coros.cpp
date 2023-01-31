/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */

#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "diag_coros.h"
#include "app.h"
#include "app_coros.h"
#include "comm_lattice.h"
#include "timer.h"
#include "error.h"
#include "memory.h"

using namespace SPPARKS_NS;
//Note: >dFE is floater, <dFe is integer
enum{inter,floater};                              // data type
enum{ZERO,FE,VACANCY,CU,NI,MN,Si,P,C,SIA};        // diagnosis terms
//enum{FE,VACANCY,CU,NI,MN,Si,P,C,SIA};
enum{hFE=11,hCU,hNI,hMN,hSi,hP,hC};               // hop steps for each element
enum{sink=30};                                    // number of sink absorption
enum{recombine=41, nmonomer, nmetalpurevac};                               // number of recombination
enum{monoFE=51, monoVACANCY, monoCU};             // number of mono-particle
enum{react=61, surffe, surfcu, bulkfe, bulkcu, nsalt, nsaltdiff};     // number of each events
enum{dFE=71,dVACANCY,dCU,dNI,dMN,dSi,dP,dC,dSIA}; // MSD for each element !! >dFE is floater
enum{energy=81,treal,fvt, metal_energy};                        // energy and realistic time
enum{cVAC=91,cFE,cVACANCY,cCU, cCE4, cCE5, cCE6, cCE7, cCE8}; 	// time averaged concentration
enum{msdFE=101,msdVACANCY,msdCU};  // MSD calculation by LC
enum{dFEx=111,dFEy, dFEz,dVACANCYx,dVACANCYy,dVACANCYz,dCUx,dCUy,dCUz};  // MSD calculation by LC
//!! be careful for the integer and float at line 164 when adding new variables
/* ---------------------------------------------------------------------- */

Diagcoros::Diagcoros(SPPARKS *spk, int narg, char **arg) : Diag(spk,narg,arg)
{
  if (strcmp(app->style,"coros") != 0)
    error->all(FLERR,"Diag_style coros requires app_style coros");

  nlist = 0;

  int iarg = iarg_child;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"list") == 0) {
      nlist = narg - iarg - 1;
      list = new char*[nlist];
      int j = 0;
      for (int i = iarg+1; i < narg; i++) {
	int n = strlen(arg[i]) + 1;
	list[j] = new char[n];
	strcpy(list[j],arg[i]);
	j++;
      }
      iarg = narg;
    } else error->all(FLERR,"Illegal diag_style coros command");
  }

  if (nlist == 0) error->all(FLERR,"Illegal diag_style coros command");
  which = new int[nlist];
  index = new int[nlist];
  itype = new int[nlist];
  ivector = new int[nlist];
  dvector = new double[nlist+1]; //nelement + total energy
}

/* ---------------------------------------------------------------------- */

Diagcoros::~Diagcoros()
{
  for (int i = 0; i < nlist; i++) delete [] list[i];
  delete [] list;
  delete [] which;
  delete [] index;
  delete [] itype;
  delete [] ivector;
  delete [] dvector;
}

/* ---------------------------------------------------------------------- */

void Diagcoros::init()
{
  appcoros = (Appcoros *) app;

  siteflag = 0;
  hopflag = 0;
  msdflag = 0;
 csiteflag = 0;


  for (int i = 0; i < nlist; i++) {
    if (strcmp(list[i],"fe") == 0) which[i] = FE; //total sites
    else if (strcmp(list[i],"vac") == 0) which[i] = VACANCY;
    else if (strcmp(list[i],"cu") == 0) which[i] = CU;
    else if (strcmp(list[i],"ni") == 0) which[i] = NI;
    else if (strcmp(list[i],"mn") == 0) which[i] = MN;
    else if (strcmp(list[i],"p") == 0) which[i] = P;
    else if (strcmp(list[i],"c") == 0) which[i] = C;
    else if (strcmp(list[i],"sia") == 0) which[i] = SIA;

    else if (strcmp(list[i],"cfe") == 0) which[i] = cFE; //time averaged concentration
    else if (strcmp(list[i],"cvac") == 0) which[i] = cVAC;
    else if (strcmp(list[i],"ccu") == 0) which[i] = cCU;
    // else if (strcmp(list[i],"cce2") == 0) which[i] = cCE2;
    // else if (strcmp(list[i],"cce3") == 0) which[i] = cCE3;
    // else if (strcmp(list[i],"cce4") == 0) which[i] = cCE4;
    // else if (strcmp(list[i],"cce5") == 0) which[i] = cCE5;
    // else if (strcmp(list[i],"cce6") == 0) which[i] = cCE6;
    // else if (strcmp(list[i],"cce7") == 0) which[i] = cCE7;
    // else if (strcmp(list[i],"cce8") == 0) which[i] = cCE8;
    else if (strcmp(list[i], "monofe") == 0) which[i] = monoFE;
    else if (strcmp(list[i], "monovac") == 0) which[i] = monoVACANCY;
    else if (strcmp(list[i], "monocu") == 0) which[i] = monoCU;
    else if (strcmp(list[i], "nmonomer") == 0) which[i] = nmonomer;
    else if (strcmp(list[i], "nmetalpurevac") == 0) which[i] = nmetalpurevac; // LC

    else if (strcmp(list[i],"hfe") == 0) which[i] = hFE;//total hop events
    else if (strcmp(list[i],"hcu") == 0) which[i] = hCU;
    else if (strcmp(list[i],"hni") == 0) which[i] = hNI;
    else if (strcmp(list[i],"hmn") == 0) which[i] = hMN;
    else if (strcmp(list[i],"hp") == 0) which[i] = hP;
    else if (strcmp(list[i],"hc") == 0) which[i] = hC;

    else if (strcmp(list[i],"dfe") == 0) which[i] = dFE;//MSD
    else if (strcmp(list[i],"dvac") == 0) which[i] = dVACANCY;
    else if (strcmp(list[i],"dcu") == 0) which[i] = dCU;
    else if (strcmp(list[i],"dni") == 0) which[i] = dNI;
    else if (strcmp(list[i],"dmn") == 0) which[i] = dMN;
    else if (strcmp(list[i],"dp") == 0) which[i] = dP;
    else if (strcmp(list[i],"dc") == 0) which[i] = dC;
    else if (strcmp(list[i],"dsia") == 0) which[i] = dSIA;

    //LC MSD calculation for each direction
    else if (strcmp(list[i],"dfex") == 0) which[i] = dFEx;
    else if (strcmp(list[i],"dfey") == 0) which[i] = dFEy;
    else if (strcmp(list[i],"dfez") == 0) which[i] = dFEz;
    else if (strcmp(list[i],"dvacx") == 0) which[i] = dVACANCYx;
    else if (strcmp(list[i],"dvacy") == 0) which[i] = dVACANCYy;
    else if (strcmp(list[i],"dvacz") == 0) which[i] = dVACANCYz;
    else if (strcmp(list[i],"dcux") == 0) which[i] = dCUx;
    else if (strcmp(list[i],"dcuy") == 0) which[i] = dCUy;
    else if (strcmp(list[i],"dcuz") == 0) which[i] = dCUz;


    //LC MSD calculation
    else if (strcmp(list[i],"msdfe") == 0) which[i] = msdFE; //MSD
    else if (strcmp(list[i],"msdvac") == 0) which[i] = msdVACANCY;
    else if (strcmp(list[i],"msdcu") == 0) which[i] = msdCU;
/*
    else if (strcmp(list[i],"mfe") == 0) which[i] = mFE;//MSD
    else if (strcmp(list[i],"mvac") == 0) which[i] = mVACANCY;
    else if (strcmp(list[i],"mcu") == 0) which[i] = mCU;
    else if (strcmp(list[i],"mni") == 0) which[i] = mNI;
    else if (strcmp(list[i],"mmn") == 0) which[i] = mMN;
    else if (strcmp(list[i],"mp") == 0) which[i] = mP;
    else if (strcmp(list[i],"mc") == 0) which[i] = mC;
    else if (strcmp(list[i],"msia") == 0) which[i] = mSIA;
*/



    // this section is to count the number of each events by LC
    else if (strcmp(list[i], "nreact") == 0) which[i] = react;
    //else if (strcmp(list[i], "nbulkfe") == 0) which[i] = bulkfe;
    //else if (strcmp(list[i], "nbulkcu") == 0) which[i] = bulkcu;
    else if (strcmp(list[i], "nsalt") == 0) which[i] = nsalt;
    else if (strcmp(list[i], "nsaltdiff") == 0) which[i] = nsaltdiff;



    else if (strcmp(list[i],"recombine") == 0) which[i] = recombine;
    else if (strcmp(list[i],"energy") == 0) which[i] = energy;
    else if (strcmp(list[i],"metal_energy") == 0) which[i] = metal_energy;
    else if (strcmp(list[i],"treal") == 0) which[i] = treal;
    else if (strcmp(list[i],"fvt") == 0) which[i] = fvt;

    else if (list[i][0] == 's' && list[i][1] == 'i' && list[i][2] == 'n' && list[i][3] == 'k') {
      int id = list[i][4] - '0';
      which[i] = sink + id;
    }

    else error->all(FLERR,"Invalid value setting in diag_style coros");
  }

  for (int i = 0; i < nlist; i++) {
    if (which[i] == FE || which[i] == VACANCY || which[i] == SIA || which[i] == CU || which[i] == NI || which[i] == MN || which[i] == P || which[i] == C)
      siteflag = 1;
    if (which[i] == hFE || which[i] == hCU || which[i] == hNI || which[i] == hMN || which[i] == hP || which[i] == hC)
      hopflag = 1;
    if(which[i] >= dFE && appcoros->diffusionflag == 1) msdflag = 1;
    if(which[i] >= dFE && which[i] <= dSIA && msdflag == 0) error->all(FLERR,"MSD calculation need displacement calculated in appcoros");
    if (which[i] >= cVAC && which[i] <= cCE8)
      csiteflag = 1;
  }

  for (int i = 0; i < nlist; i++) { ivector[i] = 0;
    dvector[i] = 0.0;
  }

  for (int i=0; i < nlist; i++) { itype[i] = inter;
    if(which[i] >= dFE) itype[i] = floater;
  }
}

/* ---------------------------------------------------------------------- */

void Diagcoros::compute()
{
  int ninter,nfloater;
  int sites[10],ivalue; // int data
  bigint nhop[10]; // big int for nhop
  int kvalue; // for number of event
  int nlocal = appcoros->nlocal;
  int nelement = appcoros->nelement;
  double dvalue ; // double data
  double *csites;
  int *monosites;
  double msd[10];
  double dir_msd[9]; // LC total MSD per direction
  double *sd; // for store sd array


  ninter = nfloater = 0;
  dvalue = 0.0;

  // time averaged concengtration
  if(csiteflag) {csites = appcoros->ct;
    //appcoros-> ct_reset(); // set ct_reset_flag = 1;
  }


  if (siteflag) {
    sites[FE] = sites[VACANCY] = sites[CU] = sites[NI] = 0;
    sites[MN] = sites[P] = sites[C] = sites[SIA] = 0;
    int *element = appcoros->element;
    for (int i = 0; i < nlocal; i++) sites[element[i]]++;
  }

  if(hopflag) {// hop event of each element
    for(int i = 1; i < nelement+1; i++) {nhop[i] = 0; nhop[i] = appcoros->hcount[i];}
  }

  if(msdflag) {// MSD calculation
    int *element = appcoros->element;
    msd[FE] = msd[VACANCY] = msd[CU] = msd[NI] = 0.0;
    msd[MN] = msd[P] = msd[C] = msd[SIA] =0.0;
    if(siteflag == 0) {//need to count total sites if not calculated above
      sites[FE] = sites[VACANCY] = sites[CU] = sites[NI] = 0;
      sites[MN] = sites[P] = sites[C] = sites[SIA] = 0;
      for (int i = 0; i < nlocal; i++) sites[element[i]]++;
    }

    // reset for dir_msd
    for (int i=0;i<9;i++){
      dir_msd[i] = 0.0;
    }

    for(int i = 0; i < nlocal; i++) {
       msd[element[i]] += appcoros->disp[3][i];

    // LC MSD calculation per direction
    //for(int i = 0; i < nlocal; i++) {
       dir_msd[3*(element[i]-1)] += appcoros->disp[0][i] * appcoros->disp[0][i];
       dir_msd[3*(element[i]-1)+1] += appcoros->disp[1][i] * appcoros->disp[1][i];
       dir_msd[3*(element[i]-1)+2] += appcoros->disp[2][i] * appcoros->disp[2][i];
    //}
    }
  }

  for (int i = 0; i < nlist; i++) {
    if (which[i] == FE) ivalue = sites[FE]; //total sites
    else if (which[i] == VACANCY) ivalue = sites[VACANCY];
    else if (which[i] == CU) ivalue = sites[CU];
    else if (which[i] == NI) ivalue = sites[NI];
    else if (which[i] == MN) ivalue = sites[MN];
    else if (which[i] == P) ivalue = sites[P];
    else if (which[i] == C) ivalue = sites[C];
    else if (which[i] == SIA) ivalue = sites[SIA];

    else if (which[i] == cFE) dvalue = csites[FE]; // time averaged concentration
    else if (which[i] == cVAC) dvalue = csites[VACANCY];
    else if (which[i] == cCU) dvalue = csites[CU];
    // else if (which[i] == cCE2) ivalue = csites[CE2];
    // else if (which[i] == cCE3) ivalue = csites[CE3];
    // else if (which[i] == cCE4) ivalue = csites[CE4];
    // else if (which[i] == cCE5) ivalue = csites[CE5];
    // else if (which[i] == cCE6) ivalue = csites[CE6];
    // else if (which[i] == cCE7) ivalue = csites[CE7];
    // else if (which[i] == cCE8) ivalue = csites[CE8];

    // count monomer for each particle from i2 by LC
    else if (which[i] == monoFE){
      appcoros-> monomer_count(); // compute monomer_count to update monomers // must have monoFe if want to count monomer
      monosites = appcoros ->monomers; // monomers array to contain monomer for each particles
      ivalue = monosites[FE]; //count vac monomer
    }
    else if (which[i] == monoVACANCY) ivalue =monosites[VACANCY]; //count vac monomer
    else if (which[i] == monoCU) ivalue = monosites[CU]; //count vac monomer
    // count monomer for pure vac from i3
    else if (which[i] == nmonomer) ivalue = appcoros ->vac_monomer_count(); //count pure vac monomer
    // count number of pure vac from i3 - PD
    else if (which[i] == nmetalpurevac) ivalue = appcoros ->metal_pure_vac_approxi(); //count pure vac monomer
    /*
    else if (which[i] == mVACANCY) ivalue = monomer_local[VACANCY];
    else if (which[i] == mCU) ivalue = monomer_local[CU];
    else if (which[i] == mNI) ivalue = monomer_local[NI];
    else if (which[i] == mMN) ivalue = monomer_local[MN];
    else if (which[i] == mP) ivalue = monomer_local[P];
    else if (which[i] == mC) ivalue = monomer_local[C];
    else if (which[i] == mSIA) ivalue = monomer_local[SIA];
    */
    else if (which[i] == hFE) ivalue = nhop[FE]; //total hop events
    else if (which[i] == hCU) ivalue = nhop[CU];
    else if (which[i] == hNI) ivalue = nhop[NI];
    else if (which[i] == hMN) ivalue = nhop[MN];
    else if (which[i] == hP) ivalue = nhop[P];
    else if (which[i] == hC) ivalue = nhop[C];

    else if (which[i] == dFE && sites[FE] > 0) dvalue = msd[FE]/sites[FE];//MSD
    else if (which[i] == dVACANCY && sites[VACANCY] > 0) dvalue = msd[VACANCY]/sites[VACANCY];
    else if (which[i] == dCU && sites[CU] > 0) dvalue = msd[CU]/sites[CU];
    else if (which[i] == dNI && sites[NI] > 0) dvalue = msd[NI]/sites[NI];
    else if (which[i] == dMN && sites[MN] > 0) dvalue = msd[MN]/sites[MN];
    else if (which[i] == dP && sites[P] > 0) dvalue = msd[P]/sites[P];
    else if (which[i] == dC && sites[C] > 0) dvalue = msd[C]/sites[C];
    else if (which[i] == dSIA && sites[SIA] > 0) dvalue = msd[SIA]/sites[SIA];

    //LC total MSD calculation for each direction
    else if (which[i] == dFEx ) dvalue = dir_msd[0]/sites[FE];//MSD
    else if (which[i] == dFEy ) dvalue = dir_msd[1]/sites[FE];
    else if (which[i] == dFEz ) dvalue = dir_msd[2]/sites[FE];
    else if (which[i] == dVACANCYx ) dvalue = dir_msd[3]/sites[VACANCY];
    else if (which[i] == dVACANCYy ) dvalue = dir_msd[4]/sites[VACANCY];
    else if (which[i] == dVACANCYz ) dvalue = dir_msd[5]/sites[VACANCY];
    else if (which[i] == dCUx ) dvalue = dir_msd[6]/sites[CU];
    else if (which[i] == dCUy ) dvalue = dir_msd[7]/sites[CU];
    else if (which[i] == dCUz ) dvalue = dir_msd[8]/sites[CU];

    // MSD calculation by LC
    else if (which[i] == msdFE){
      appcoros-> MSD_calculation(); // compute monomer_count to update monomers // must have monoFe if want to count monomer
      sd = appcoros ->sd; // monomers array to contain monomer for each particles
      dvalue = sd[FE]/sites[FE]; //count vac monomer
    }
    else if (which[i] == msdVACANCY) dvalue = sd[VACANCY]/sites[VACANCY];
    else if (which[i] == msdCU) dvalue = sd[CU]/sites[CU];

    // this section is to count the number of each events by LC
    else if (which[i] == react) ivalue = appcoros->nreact; //total reaction event
    //else if (which[i] == bulkfe) ivalue = appcoros->nbulkfe; //total bulk_diff event for id2 = 1
    //else if (which[i] == bulkcu) ivalue = appcoros->nbulkcu; //total bulk_diff event for id2 = 3
    else if (which[i] == nsalt) ivalue = appcoros->count_salt(); //total salt particle for id3 = 1
    else if (which[i] == nsaltdiff) ivalue = appcoros->num_saltdiffusion; //total bulk_diff event for id2 = 3

    else if (which[i] == recombine) dvalue = appcoros->nrecombine; //number of reocmbination
    else if (which[i] == energy) dvalue = appcoros->total_energy(); //system energy
    else if (which[i] == metal_energy) dvalue = appcoros->total_metal_energy(); //system energy
    else if (which[i] == treal) dvalue = appcoros->realtime; //realistic time
    else if (which[i] == fvt) dvalue = appcoros->fvt; //realistic time

    else if (which[i] > sink && which[i] < recombine) {
      int id = which[i] - sink;
      ivalue = appcoros->nabsorption[id-1];
    }

    if(which[i] >= dFE) nfloater++;
    else ninter++;

    MPI_Allreduce(&ivalue,&ivector[ninter],1,MPI_INT,MPI_SUM,world);
    MPI_Allreduce(&dvalue,&dvector[nfloater],1,MPI_DOUBLE,MPI_SUM,world);
  }
}

/* ---------------------------------------------------------------------- */

void Diagcoros::stats(char *str)
{
  int ninter, nfloater;

  ninter = nfloater = 0;
  for (int i = 0; i < nlist; i++) {
    if(itype[i] == inter) { ninter++;
      sprintf(str," %d",ivector[ninter]);
    }
    else if(itype[i] == floater) { nfloater++;
      sprintf(str," %14.8g",dvector[nfloater]);
    }
    str += strlen(str);
  }
  appcoros->KMC_stop(); // by LC
}

/* ---------------------------------------------------------------------- */

void Diagcoros::stats_header(char *str)
{
  for (int i = 0; i < nlist; i++) {
    sprintf(str," %s",list[i]);
    str += strlen(str);
  }
}
