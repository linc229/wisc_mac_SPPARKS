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
#include "diag_seg.h"
#include "app.h"
#include "app_seg.h"
#include "comm_lattice.h"
#include "timer.h"
#include "error.h"
#include "memory.h"

using namespace SPPARKS_NS;

enum{inter,floater};                              		// data type
enum{VAC=0,INT,CE1,CE2,CE3,CE4,CE5,CE6,CE7,CE8};   		// CE: chemical element
enum{hVAC=10,hINT,hCE1,hCE2,hCE3,hCE4,hCE5,hCE6,hCE7,hCE8}; 	// hop steps for each element
enum{sink=20};                                    		// number of sink absorption
enum{recombine=30,FPair,selfion,vabsorb,iabsorb};      		// number of recombination
enum{mono=40,di_mono1,di_mono2,di_mono3,di_mono4,di_mono5,di_mono6,di_mono7,di_mono8}; // Peng: number of dissolved particles
enum{sia=50};                                     		// onsager coefficient, floats start, integers should be set before this line
enum{cVAC=60,cINT,cCE1,cCE2,cCE3,cCE4,cCE5,cCE6,cCE7,cCE8}; 	// time averaged concentration
enum{dVAC=70,dINT,dCE1,dCE2,dCE3,dCE4,dCE5,dCE6,dCE7,dCE8}; 	// MSD
enum{energy=80,treal,fvt};                        		// energy and realistic time
enum{ris=90,};                                     		// number of ris
enum{lij=100};                                     		// onsager coefficient
enum{sro=200};                                     		// short range order
/* ---------------------------------------------------------------------- */

DiagSeg::DiagSeg(SPPARKS *spk, int narg, char **arg) : Diag(spk,narg,arg)
{
  if (strcmp(app->style,"seg") != 0)
    error->all(FLERR,"Diag_style seg requires app_style seg");

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
    } else error->all(FLERR,"Illegal diag_style seg command");
  }

  if (nlist == 0) error->all(FLERR,"Illegal diag_style seg command");
  which = new int[nlist];
  index = new int[nlist];
  itype = new int[nlist];
  ivector = new int[nlist];
  dvector = new double[nlist+1]; //nelement + total energy
}

/* ---------------------------------------------------------------------- */

DiagSeg::~DiagSeg()
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

void DiagSeg::init()
{
  appseg = (AppSeg *) app;

  siteflag = 0;
  csiteflag = 0;
  hopflag = 0;
  msdflag = 0;
  siaflag = 0;
  risflag = appseg->seg_flag;;

  for (int i = 0; i < nlist; i++) {
    if (strcmp(list[i],"vac") == 0) which[i] = VAC; //total sites
    else if (strcmp(list[i],"int") == 0) which[i] = INT;
    else if (strcmp(list[i],"ce1") == 0) which[i] = CE1;
    else if (strcmp(list[i],"ce2") == 0) which[i] = CE2;
    else if (strcmp(list[i],"ce3") == 0) which[i] = CE3;
    else if (strcmp(list[i],"ce4") == 0) which[i] = CE4;
    else if (strcmp(list[i],"ce5") == 0) which[i] = CE5;
    else if (strcmp(list[i],"ce6") == 0) which[i] = CE6;
    else if (strcmp(list[i],"ce7") == 0) which[i] = CE7;
    else if (strcmp(list[i],"ce8") == 0) which[i] = CE8;

    else if (strcmp(list[i],"cvac") == 0) which[i] = cVAC; //time averaged concentration
    else if (strcmp(list[i],"cint") == 0) which[i] = cINT;
    else if (strcmp(list[i],"cce1") == 0) which[i] = cCE1;
    else if (strcmp(list[i],"cce2") == 0) which[i] = cCE2;
    else if (strcmp(list[i],"cce3") == 0) which[i] = cCE3;
    else if (strcmp(list[i],"cce4") == 0) which[i] = cCE4;
    else if (strcmp(list[i],"cce5") == 0) which[i] = cCE5;
    else if (strcmp(list[i],"cce6") == 0) which[i] = cCE6;
    else if (strcmp(list[i],"cce7") == 0) which[i] = cCE7;
    else if (strcmp(list[i],"cce8") == 0) which[i] = cCE8;

   // hopping calculations (Peng)
    else if (strcmp(list[i],"hvac") == 0) which[i] = hVAC; //total hopping event n
    else if (strcmp(list[i],"hint") == 0) which[i] = hINT;
    else if (strcmp(list[i],"hce1") == 0) which[i] = hCE1;
    else if (strcmp(list[i],"hce2") == 0) which[i] = hCE2;
    else if (strcmp(list[i],"hce3") == 0) which[i] = hCE3;
    else if (strcmp(list[i],"hce4") == 0) which[i] = hCE4;
    else if (strcmp(list[i],"hce5") == 0) which[i] = hCE5;
    else if (strcmp(list[i],"hce6") == 0) which[i] = hCE6;
    else if (strcmp(list[i],"hce7") == 0) which[i] = hCE7;
    else if (strcmp(list[i],"hce8") == 0) which[i] = hCE8;


    else if (strcmp(list[i],"dvac") == 0) which[i] = dVAC; //MSD
    else if (strcmp(list[i],"dint") == 0) which[i] = dINT;
    else if (strcmp(list[i],"dce1") == 0) which[i] = dCE1;
    else if (strcmp(list[i],"dce2") == 0) which[i] = dCE2;
    else if (strcmp(list[i],"dce3") == 0) which[i] = dCE3;
    else if (strcmp(list[i],"dce4") == 0) which[i] = dCE4;
    else if (strcmp(list[i],"dce5") == 0) which[i] = dCE5;
    else if (strcmp(list[i],"dce6") == 0) which[i] = dCE6;
    else if (strcmp(list[i],"dce7") == 0) which[i] = dCE7;
    else if (strcmp(list[i],"dce8") == 0) which[i] = dCE8;

   /*// time correction due to trapping disable currently
    else if (strcmp(list[i],"mvac") == 0) which[i] = mVAC; /mononmer
    else if (strcmp(list[i],"mint") == 0) which[i] = mINT;
    else if (strcmp(list[i],"mce1") == 0) which[i] = mCE1;
    else if (strcmp(list[i],"mce2") == 0) which[i] = mCE2;
    else if (strcmp(list[i],"mce3") == 0) which[i] = mCE3;
    else if (strcmp(list[i],"mce4") == 0) which[i] = mCE4;
    else if (strcmp(list[i],"mce5") == 0) which[i] = mCE5;
    else if (strcmp(list[i],"mce6") == 0) which[i] = mCE6;
    else if (strcmp(list[i],"mce7") == 0) which[i] = mCE7;
    else if (strcmp(list[i],"mce8") == 0) which[i] = mCE8;

    else if (strcmp(list[i],"treal") == 0) which[i] = treal;
    else if (strcmp(list[i],"fvt") == 0) which[i] = fvt;
   */

   //Peng: dissolved B particles, monomers and dimers
    else if (strcmp(list[i],"di_mono1") == 0) which[i] = di_mono1;
    else if (strcmp(list[i],"di_mono2") == 0) which[i] = di_mono2;
    else if (strcmp(list[i],"di_mono3") == 0) which[i] = di_mono3;
    else if (strcmp(list[i],"di_mono4") == 0) which[i] = di_mono4;
    else if (strcmp(list[i],"di_mono5") == 0) which[i] = di_mono5;
    else if (strcmp(list[i],"di_mono6") == 0) which[i] = di_mono6;
    else if (strcmp(list[i],"di_mono7") == 0) which[i] = di_mono7;
    else if (strcmp(list[i],"di_mono8") == 0) which[i] = di_mono8;

    else if (strcmp(list[i],"recombine") == 0) which[i] = recombine;
    else if (strcmp(list[i],"nfp") == 0) which[i] = FPair;
    else if (strcmp(list[i],"selfion") == 0) which[i] = selfion;
    else if (strcmp(list[i],"vabsorb") == 0) which[i] = vabsorb;
    else if (strcmp(list[i],"iabsorb") == 0) which[i] = iabsorb;
    else if (strcmp(list[i],"energy") == 0) which[i] = energy;

    // number of each chemical type of sias
    else if (list[i][0] == 's' && list[i][1] == 'i' && list[i][2] == 'a') {
      int id = list[i][3] - '0';
      which[i] = sia + id;
      siaflag = 1;
    }

    // sink absorption
    else if (list[i][0] == 's' && list[i][1] == 'i' && list[i][2] == 'n' && list[i][3] == 'k') {
      int id = list[i][4] - '0';
      which[i] = sink + id;
    }

    // ris of element i
    else if (list[i][0] == 'r' && list[i][1] == 'i' && list[i][2] == 's') {
      int id = list[i][3] - '0';
      which[i] = ris + id;
    }

    // onsager coefficients
    else if (list[i][0] == 'l' && list[i][1] == 'i' && list[i][2] == 'j') {
      int id1 = list[i][3] - '0';
      int id2 = list[i][4] - '0';
      which[i] = lij + id1*10 + id2;
    }

    // short range order
    else if (list[i][0] == 's' && list[i][1] == 'r' && list[i][2] == 'o') {
      int id1 = list[i][3] - '0';
      int id2 = list[i][4] - '0';
      which[i] = sro + id1*10 + id2;
    }

    else error->all(FLERR,"Invalid value setting in diag_style seg");
  }

  for (int i = 0; i < nlist; i++) {
    if (which[i] >= VAC && which[i] <= CE8)
      siteflag = 1;
    if (which[i] >= cVAC && which[i] <= cCE8)
      csiteflag = 1;
    if (which[i] >= dVAC && which[i] <= dCE8 && appseg->diffusionflag == 1)
       msdflag = 1;
  }

  for (int i = 0; i < nlist; i++) {
    ivector[i] = 0;
    dvector[i] = 0.0;
  }

  for (int i=0; i < nlist; i++) { itype[i] = inter;
    if(which[i] >= sia) itype[i] = floater;
  }
}

/* ---------------------------------------------------------------------- */

void DiagSeg::compute()
{
  int i,ninter,nfloater;
  int sites[10],nhop[10],ivalue; // int data
  int nlocal = appseg->nlocal;
  int nelement = appseg->nelement;
  int m; //count_dumbbell
  double dvalue; // double data
  double *csites;
  double msd[10];

  ninter = nfloater = 0;
  dvalue = 0.0;

  // time dependent segregation profile
  if (risflag) {appseg->ris_time();}
  // time averaged concengtration
  if(csiteflag) {csites = appseg->ct;}
  // number of each chemical type of sias
  if(siaflag) {appseg->count_dumbbell(m);}


  // site information
  if (siteflag || msdflag) {
    for(i=VAC; i<nelement; i++) sites[i] = 0;
    int *element = appseg->element;
    for(i = 0; i < nlocal; i++) sites[element[i]]++;
  }


  if(hopflag) {// hop event of each element (Peng)
    for(int i = 1; i < nelement+1; i++) {nhop[i] = 0; nhop[i] = appseg->hcount[i];}
  }


  if(msdflag) {// MSD calculation
    int *element = appseg->element;
    for(i=VAC; i<nelement; i++) msd[i] = 0.0;

    for(int i = 0; i < nlocal; i++) {
       msd[element[i]] += appseg->disp[3][i];
    }
  }

  for (int i = 0; i < nlist; i++) {
    if (which[i] == VAC) ivalue = sites[VAC]; //total sites
    else if (which[i] == INT) ivalue = sites[INT];
    else if (which[i] == CE1) ivalue = sites[CE1];
    else if (which[i] == CE2) ivalue = sites[CE2];
    else if (which[i] == CE3) ivalue = sites[CE3];
    else if (which[i] == CE4) ivalue = sites[CE4];
    else if (which[i] == CE5) ivalue = sites[CE5];
    else if (which[i] == CE6) ivalue = sites[CE6];
    else if (which[i] == CE7) ivalue = sites[CE7];
    else if (which[i] == CE8) ivalue = sites[CE8];

    else if (which[i] == cVAC) ivalue = csites[VAC]; // time averaged concentration
    else if (which[i] == cINT) ivalue = csites[INT];
    else if (which[i] == cCE1) ivalue = csites[CE1];
    else if (which[i] == cCE2) ivalue = csites[CE2];
    else if (which[i] == cCE3) ivalue = csites[CE3];
    else if (which[i] == cCE4) ivalue = csites[CE4];
    else if (which[i] == cCE5) ivalue = csites[CE5];
    else if (which[i] == cCE6) ivalue = csites[CE6];
    else if (which[i] == cCE7) ivalue = csites[CE7];
    else if (which[i] == cCE8) ivalue = csites[CE8];

    else if (which[i] == dVAC && sites[VAC] > 0) dvalue = msd[VAC]/sites[VAC]; // MSD
    else if (which[i] == dINT && sites[INT] > 0) dvalue = msd[INT]/sites[INT];
    else if (which[i] == dCE1 && sites[CE1] > 0) dvalue = msd[CE1]/sites[CE1];
    else if (which[i] == dCE2 && sites[CE2] > 0) dvalue = msd[CE2]/sites[CE2];
    else if (which[i] == dCE3 && sites[CE3] > 0) dvalue = msd[CE3]/sites[CE3];
    else if (which[i] == dCE4 && sites[CE4] > 0) dvalue = msd[CE4]/sites[CE4];
    else if (which[i] == dCE5 && sites[CE5] > 0) dvalue = msd[CE5]/sites[CE5];
    else if (which[i] == dCE6 && sites[CE6] > 0) dvalue = msd[CE6]/sites[CE6];
    else if (which[i] == dCE7 && sites[CE7] > 0) dvalue = msd[CE7]/sites[CE7];
    else if (which[i] == dCE8 && sites[CE8] > 0) dvalue = msd[CE8]/sites[CE8];

    else if (which[i] == energy) dvalue = appseg->total_energy(); //system energy
    else if (which[i] > sink && which[i] < recombine)  ivalue = 0; //to be added later
    else if (which[i] == FPair) ivalue = appseg->nFPair; //number of reocmbination
    else if (which[i] == selfion) ivalue = appseg->nself_ion; //number of selfionn
    else if (which[i] == recombine) ivalue = appseg->nrecombine[VAC]; //number of reocmbination

//Peng: number of solute particles (ce1 as solvent)
    else if (which[i] == di_mono1) {// ce2
      int n = 3;
      ivalue = appseg->solubility(n);
    }
    else if (which[i] == di_mono2) {// ce3
      int n = 4;
      ivalue = appseg->solubility(n);
    }
    else if (which[i] == di_mono3) {// ce4
      int n = 5;
      ivalue = appseg->solubility(n);
    }
    else if (which[i] == di_mono4) {// ce5
      int n = 6;
      ivalue = appseg->solubility(n);
    }
    else if (which[i] == di_mono5) {// ce6
      int n = 7;
      ivalue = appseg->solubility(n);
    }
    else if (which[i] == di_mono6) {// ce7
      int n = 8;
      ivalue = appseg->solubility(n);
    }
    else if (which[i] == di_mono7) {// ce8
      int n = 9;
      ivalue = appseg->solubility(n);
    }


    else if (which[i] == vabsorb) { // absorbed vacan
      int nabsorb = 0;
      int nsink = appseg->nsink;
      for(int i = 1; i < nsink+1; i++) { // sink id starts from 1
         nabsorb += appseg->nabsorption[VAC][i];
      }
      ivalue = nabsorb;
    }

    else if (which[i] == iabsorb) { // absorbed sia
      int nabsorb = 0;
      int nsink = appseg->nsink;
      for(int i = 1; i < nsink+1; i++) {
         nabsorb += appseg->nabsorption[INT][i];
      }
      ivalue = nabsorb;
    }

    else if (which[i] >= sia && which[i] < cVAC) { // # of sias
      int id = which[i] - sia;
      dvalue = appseg->nsia[id];
    }

    else if (which[i] >= ris && which[i] < lij) { // ris
      int id = which[i] - ris;
      dvalue = appseg->ris_total[id];
    }

    else if (which[i] >= lij && which[i] < sro) { // onsager
      int id2 = (which[i] - lij)%10;
      int id1 = (which[i] - lij - id2)/10;
      dvalue = appseg->Lij[id1][id2]; // calcualted in appseg->onsager()
    }
    else if (which[i] >= sro) { // short range order
      int id2 = (which[i] - sro)%10;
      int id1 = (which[i] - sro - id2)/10;
      appseg->short_range_order(); // compute short range order
      dvalue = appseg->sro[id1][id2];
    }

    if(which[i] >= sia) nfloater++;
    else ninter++;
    MPI_Allreduce(&ivalue,&ivector[ninter],1,MPI_INT,MPI_SUM,world);
    MPI_Allreduce(&dvalue,&dvector[nfloater],1,MPI_DOUBLE,MPI_SUM,world);
  }
}

/* ---------------------------------------------------------------------- */

void DiagSeg::stats(char *str)
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
}

/* ---------------------------------------------------------------------- */

void DiagSeg::stats_header(char *str)
{
  for (int i = 0; i < nlist; i++) {
    sprintf(str," %s",list[i]);
    str += strlen(str);
  }
}
