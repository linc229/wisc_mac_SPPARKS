/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
   See the README file in the top-level SPPARKS directory.
*************************************************************************************
   This application does the vacancy diffusion includes surface diffusion
   , and impurity diffusion in molten salt via hop in a FCC lattice in NiCr alloy system
   Barriers for diffusion was calcualted by FISE model.
   bond energies. Contributer: Yongfeng Zhang, yongfeng.zhang@inl.gov
------------------------------------------------------------------------- */

#include "math.h"
#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "app_coros.h"
#include "domain.h"
#include "solve.h"
#include "random_park.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"

#include "finish.h" // added for KMC_stop by LC
#include <iostream>
#include <fstream>
using namespace std;
#include "app_lattice.h"

using namespace SPPARKS_NS;

enum{NOOP,BCC,NBCC,SALT};                          // all sites are BCC except for SIAs; type
enum{ZERO,FE,VACANCY,CU,NI,MN,Si,P,C,SIA};       // same as Diagcoros; element

#define DELTAEVENT 100000
#define MAX2NN 6 // max 2NN of BCC lattice
#define BIGNUMBER 1e18 // define a big number

/* ---------------------------------------------------------------------- */

Appcoros::Appcoros(SPPARKS *spk, int narg, char **arg) :
  AppLattice(spk,narg,arg)
{
  ninteger = 3; // first for lattice type,second for element
  ndouble = 0; // default
  //ndouble = 1; // add 1 to contain potential by LC
  delpropensity = 2; // LC note: control layer for nborder ref value: 2
  delevent = 1;
  allow_kmc = 1;
  allow_rejection = 0;

  engstyle = 1; //1 for 1NN interaction, 2 for 2NN interaction; default 1
  diffusionflag = 0; //flag for MSD calculations, 1 for yes, 0 for no; default 0
  //seg_flag = 0; // LC
  eisink_flag = 0; // LC

  if (narg < 1) error->all(FLERR,"Illegal app_style command");
  if (narg >= 2) engstyle = atoi(arg[1]);
  if (narg >= 3) diffusionflag = atoi(arg[2]);
  if (narg >= 4) concentrationflag = atoi(arg[3]);
  if (engstyle == 2) delpropensity += 1;// increase delpropensity for 2NN interaction

  // darray 1-4 for msd if activated, followed by concentrations
  if (diffusionflag == 1) {
    ninteger++;
    ndouble += 4;
  }
  // calculate concentration fiels for certain elements
  if (concentrationflag) {ndouble += concentrationflag + 1;
  ndiff = diffusionflag*4;}

  create_arrays();

  firsttime = 1;
  esites = NULL;
  echeck = NULL;
  events = NULL;
  maxevent = 0;
  firstevent = NULL;

  // set random number generator
  rancoros = new RandomPark(ranmaster->uniform());
  double seed = ranmaster->uniform();
  rancoros->reset(seed,me,100);

  // flags for bond interactions
  ebond1 = NULL;
  ebond2 = NULL;
  mbarrier = NULL;
  hcount = NULL; //numner of vacancy switch events
  nn1flag = nn2flag = barrierflag = time_flag = 0; //flags for bond energy and migration barriers
  // flag and parameter for attempt frequency by LC
  attemptfrequencyflag = 0;
  attemptfrequency = NULL;

  // reset hcount by LC
  //for (int i=0; i<=nelement; i++){
  //  hcount[i] = 0;
  //}


  // flags and parameters for sinks, dislocations, reactions and ballistic mixing
  sink_flag = elastic_flag = moduli_flag = dislocation_flag = reaction_flag = acceleration_flag = 0; //flags for sink dislocation and vacancy
  nsink = ndislocation = nreaction = nballistic = ntrap = 0;
  nsaltdiffusion = 0; // LC

  // arrays for dislocations
  dislocation_type = line_vector = nsegment = NULL;
  burgers = xdislocation = NULL;
  dislocation_radius = NULL;

  // arrays for sinks
  // sink_shape = sink_type = sink_normal = sink_segment = nabsorption = NULL;
  // sink_strength = sink_radius = sink_mfp = NULL;
  // xsink = NULL;
  isink = sink_shape = sink_normal = sink_segment = NULL;
  sink_range = sink_radius = ci = sink_dr = sink_dt = sink_dt_new = sink_dt_old =  NULL;
  xsink = sink_mfp = NULL;
  sinksite = NULL;
  nabsorption = nreserve = sinkid = NULL;

  // arrays for reactions
  rsite = rinput = routput = rcount = renable = rtarget = NULL;
  nsites_local = target_local = target_global =NULL; //number of reactions
  rbarrier = rrate = NULL;

  // arrays for ballistic mixing
  rdamp = pn_local = pn_global = NULL;
  bfreq = NULL;
  time_old = time_new = NULL;
  xmix = pmix = NULL;
  min_bfreq = BIGNUMBER;

  // arrays for salt diffusion by LC
  salt_bfreq = NULL;
  salt_time_old = salt_time_new = NULL;
  num_saltdiffusion = 0;
  temp_salt_time = 1;

  // 2NN neigbor information
  numneigh2 = NULL;
  neighbor2 = NULL;

  // number of recombinations
  nrecombine = 0;

  // number of each events by LC
  nreact = 0;
  //nbulkfe = 0; // comment but keep it
  //nbulkcu = 0;

  // parameter for barrier_extract  by LC
  extract_flag = 0;
  evap_extract_flag = 0;
  np_extract_flag = 0;
  //ct_reset_flag = 0;
  //ct_site_flag = 0;
  large_propensity_extract_flag = 0;
  //dump_event_flag = 0;

  // parameter for total_record by LC
  total_Ni = 0;
  total_vac = 0;
  total_Cr = 0;
  surface_Ni = 0;
  surface_Cr = 0;

  // parameter for dump_event by LC
  dump_index=0;
  dump_event_all_index = 0;
  dump_event_all_file = 0;
  dump_event_all_flag = 0;


  // corrosion_info_output parameter by LC
  coros_index = 0;
  coros_time = NULL;
  coros_id = NULL;
}

/* ---------------------------------------------------------------------- */

Appcoros::~Appcoros()
{
  delete [] esites;
  delete [] echeck;
  delete [] hcount; // number of vacancy switch
  delete rancoros;

  memory->sfree(events);
  memory->destroy(firstevent);
  memory->destroy(ebond1);
  memory->destroy(ebond2);
  memory->destroy(mbarrier);
  memory->destroy(attemptfrequency); // attempt frequency by LC
  memory->destroy(nsites_local);

  if (engstyle == 2) {// memory use for 2NNs
    memory->destroy(numneigh2);
    memory->destroy(neighbor2);
  }

  if (dislocation_flag) {// memory use related to dislocation
    memory->destroy(stress);
    memory->destroy(dislocation_type);
    memory->destroy(burgers);
    memory->destroy(dislocation_radius);
    memory->destroy(line_vector);
    memory->destroy(nsegment);
    memory->destroy(xdislocation);
  }

  if (sink_flag) { // memory use related to sink
    // memory->destroy(sink_shape);
    // memory->destroy(sink_type);
    // memory->destroy(sink_strength);
    // memory->destroy(xsink);
    // memory->destroy(isink);
    // memory->destroy(sink_normal);
    // memory->destroy(sink_segment);
    // memory->destroy(sink_radius);
    // memory->destroy(sink_mfp);
    // memory->destroy(nabsorption);
    memory->destroy(sink_shape);
    memory->destroy(sink_range);
    memory->destroy(xsink);
    memory->destroy(isink);
    memory->destroy(sinksite);
    memory->destroy(sinkid);
    memory->destroy(sink_normal);
    memory->destroy(sink_segment);
    memory->destroy(sink_radius);
    memory->destroy(eisink);
    memory->destroy(sink_dt_new);
    memory->destroy(sink_dt_old);
    memory->destroy(sink_dr);
    memory->destroy(nabsorption);
    memory->destroy(nreserve);
  }

  if (ballistic_flag) { // memory use related to ballistic mixing
    memory->destroy(bfreq);
    memory->destroy(time_old);
    memory->destroy(time_new);
    memory->destroy(rdamp);
    memory->destroy(pn_local);
    memory->destroy(pn_global);
    memory->destroy(xmix);
    memory->destroy(pmix);
  }

  // salt diffusion by LC
  if (saltdiffusion_flag) { // memory use related to salt diffusion
    memory->destroy(salt_bfreq);
    memory->destroy(salt_time_old);
    memory->destroy(salt_time_new);
  }

  if (reaction_flag) {// memory use related to reaction
    memory->destroy(rsite);
    memory->destroy(rinput);
    memory->destroy(routput);
    memory->destroy(rcount);
    memory->destroy(renable);
    memory->destroy(rtarget);
    memory->destroy(rbarrier);
    memory->destroy(rrate);
    memory->destroy(target_local);
    memory->destroy(target_global);
  }

  if (acceleration_flag) {// memory use for acceleraton
    memory->destroy(trap_type);
  }

  if (coros_flag) {  // corrosion_info_output
    memory->destroy(coros_id);
    memory->destroy(coros_time);
  }

  // bulk_vac memory by lC
  memory->destroy(zprofile);

}

/* ----------------------------------------------------------------------
  setup bond energies for 1NN, 2NN, and SPs
------------------------------------------------------------------------- */

void Appcoros::input_app(char *command, int narg, char **arg)
{
  int i,j,ibond;
  int nlattice = nlocal + nghost;

  // 1NN bond energy taken in the order of 11 12 ... 1N; 22 ... 2N; ...; NN
  if (strcmp(command,"ebond1") == 0) {

    if (narg < 3) error->all(FLERR,"Illegal ebond1 command");
    nelement = atoi(arg[0]);   // num of elements

    memory->create(nsites_local,nelement,"app/coros:nsites_local");
    memory->create(ebond1,nelement+1,nelement+1,"app/coros:ebond1");
    memory->create(ct,nelement+1,"app/coros:ct"); //time averaged concentration based on the fractional occupation at each site
    memory->create(ct_new,nelement+1,"app/coros:ct_new"); //time averaged concentration
    memory->create(monomers,nelement,"app/coros:monomers");
    memory->create(sd,nelement,"app/coros:sd"); // sd create for msd calculation
    //if(concentrationflag) memory->create(ct_site,nelement,nlocal,"app/coros:ct_site"); //site coefficient
    memory->create(ct_site,nelement+1,nlocal,"app/coros:ct_site"); //LC
    memory->create(ct_site_new,nelement+1,nlocal,"app/coros:ct_site_new");  //LC

    memory->create(ct_site_temp_t,nlocal,"app/coros:ct_site_temp_t");  //LC
    memory->create(ct_site_temp_i2,nlocal,"app/coros:ct_site_temp_i2");  //LC

    memory->create(site_time,nlocal,"app/coros:site_time");  //LC
    memory->create(site_time_interval,8,"app/coros:site_time_interval");  //LC

    memory->create(i3_site,5,nlocal,"app/coros:i3_site"); // LC
    memory->create(i3_site_new,5,nlocal,"app/coros:i3_site_new"); // LC

    // diag version array
    memory->create(ct_site_diag,nelement+1,nlocal,"app/coros:ct_site_diag"); //LC
    memory->create(ct_site_new_diag,nelement+1,nlocal,"app/coros:ct_site_new_diag");  //LC
    memory->create(ct_site_temp_t_diag,nlocal,"app/coros:ct_site_temp_t_diag");  //LC
    memory->create(ct_site_temp_i2_diag,nlocal,"app/coros:ct_site_temp_i2_diag");  //LC
    memory->create(site_time_diag,nlocal,"app/coros:site_time_diag");  //LC

    hcount = new double [nelement+1]; // total numner of switching with a vacancy;

    if(narg != nelement*(nelement+1)/2+1) error->all(FLERR,"Illegal ebond command");
    nn1flag = 1;

    for (i = 1; i <= nelement; i++ ) {
      for (j = i; j <= nelement; j++ ) {
        ibond = ibonde(i,j,nelement);
        ebond1[i][j] = atof(arg[ibond]);
        if(j > i) ebond1[j][i] = ebond1[i][j];
      }
    }
  }

  // 2NN bond energy taken in the order of 11 12 ... 1N; 22 ... 2N; ...; NN
  else if (strcmp(command,"ebond2") == 0) {

    nelement = atoi(arg[0]);   // num of elements
    memory->create(ebond2,nelement+1,nelement+1,"app/coros:ebond2");
    if(narg != nelement*(nelement+1)/2+1) error->all(FLERR,"Illegal ebond command");

    nn2flag = 1;
    for (i = 1; i <= nelement; i++ ) {
      for (j = i; j <= nelement; j++ ) {
        ibond = ibonde(i,j,nelement);
        ebond2[i][j] = atof(arg[ibond]);
        if (j > i) ebond2[j][i] = ebond2[i][j];
      }
    }
  }

  // migration barriers for each element
  else if (strcmp(command, "migbarrier") ==0) {
    if (narg < 2 || narg % 2 != 0) error->all(FLERR,"Illegal migbarrier command");
    barrierflag = 1;
    memory->create(mbarrier,nelement+1,"app/coros:mbarrier");
    for (i=0; i<narg-1; i++) {
      if(i % 2 == 0){ j = atoi(arg[i]);
        mbarrier[j] = atof(arg[i+1]);
      }
    }
  }

// attempt frequency for vacancy diffusion for each particle by LC
// the frequency can also used for time scale conversion...
// the unit is Thz = 10^12 Hz = 10^12 1/s = 1/ps
  else if (strcmp(command, "attemptfrequency") ==0) {
  if (narg < 2 || narg % 2 != 0) error->all(FLERR,"Illegal attempt frequency command");
  attemptfrequencyflag = 1;
  memory->create(attemptfrequency,nelement+1,"app/coros:attemptfrequency");
  for (i=0; i<narg-1; i++) {
    if(i % 2 == 0){ j = atoi(arg[i]);
      attemptfrequency[j] = atof(arg[i+1]);
    }
  }
}

  // time intervals used to estimate solute trapping
  /*
  else if (strcmp(command, "time_tracer") ==0) {
    if (narg < 1 ) error->all(FLERR,"Illegal time_tracer command");
    time_flag = 1;
    dt_interval = atof(arg[0]);
  }
  // type of solutes that trap vacancy and need acceleration
  else if (strcmp(command, "acceleration") ==0) {
    acceleration_flag = 1;
    memory->create(trap_type,nelement+1,"app/coros:trap_type");
    if (narg < 1 ) error->all(FLERR,"Illegal acceleration command");
    ntrap = 0;
    while(ntrap < narg) {
      trap_type[ntrap] = atoi(arg[ntrap]);
      ntrap ++;
    }
  }
  // elastic moduli for stress calculation
  else if (strcmp(command, "moduli") ==0) {
    moduli_flag = 1;
    if(narg != 5) error->all(FLERR,"illegal moduli command");
    c11 = atof(arg[0]);
    c12 = atof(arg[1]);
    c44 = atof(arg[2]);
    ninteg = atoi(arg[3]);
    dcore = atof(arg[4]);
  }
  // calculating stress field of dislocations or loops
  else if (strcmp(command, "dislocation") ==0) {
    if(narg < 8 || moduli_flag == 0) error->all(FLERR,"illegal dislocation command");
    if(dislocation_flag == 0) {
      memory->create(stress,nlattice,6,"app/coros:stress");
      for(i = 0; i < nlattice; i++) {
        for(j = 0; j < 6; j++) stress[i][j] = 0;
      }
    }
    dislocation_flag = 1;
    grow_dislocations();
    // dislocation type, 1 straight, others loop
    dislocation_type[ndislocation] = atoi(arg[0]);
    burgers[ndislocation][0] = atof(arg[1]); //burgers vector
    burgers[ndislocation][1] = atof(arg[2]);
    burgers[ndislocation][2] = atof(arg[3]);
    xdislocation[ndislocation][0] = atof(arg[4]); //center of location or loop
    xdislocation[ndislocation][1] = atof(arg[5]);
    xdislocation[ndislocation][2] = atof(arg[6]);
    line_vector[ndislocation] = atoi(arg[7]); //line vector or plane normal for loop, X(0) or Y(1) or Z(2) only
    // extra parameters for loops
    if(narg > 8) {
      dislocation_radius[ndislocation] = atof(arg[8]); //loop radius, for loop only
      nsegment[ndislocation] = atoi(arg[9]); //loop segments(shape), for loop only
    }
    // compute the stress field
    stress_field(ndislocation);
    ndislocation ++;
  }
  // elastic interaction with stress field
  else if (strcmp(command, "elastic_interaction") ==0) {
    elastic_flag = 1;
    if(narg < 2 || dislocation_flag == 0) error->all(FLERR,"illegal elastic interaction command");
    int iarg = 0;
    int itype = 0;
    for (i = 1; i <= nelement; i++ ) evol[i] = 0.0;
    while(iarg < narg) {
      itype = atoi(arg[iarg]);
      evol[itype] = atof(arg[iarg+1]);
      iarg += 2;
    }
  }*/
  // define sinks to defects, could be dislocations or interfaces or 3D regions
  // else if (strcmp(command, "sink") ==0) {
  //   if(narg != 10) error->all(FLERR,"illegal sink command");
  //   if(sink_flag == 0) {
  //     memory->create(isink, nlattice, nelement,"app/coros:isink");
  //     for(i = 0; i < nlattice; i++) {
  //       for(j = 0; j < nelement; j++)  isink[i][j] = 0;
  //     }
  //   }
  //   sink_flag = 1;
  //   grow_sinks();
  //   sink_type[nsink] = atoi(arg[0]); // sink to certain element
  //   sink_strength[nsink] = atof(arg[1]); // thickness of sink
  //   sink_shape[nsink] = atoi(arg[2]); // 1 dislocation, 2 interface, 3 3D region
  //   xsink[nsink][0] = atof(arg[3]); // coordinaiton of sink center
  //   xsink[nsink][1] = atof(arg[4]);
  //   xsink[nsink][2] = atof(arg[5]);
  //   sink_normal[nsink] = atoi(arg[6]); // normal of planar sinks
  //   sink_radius[nsink] = atof(arg[7]); // radius for circular or polygonal sinks
  //   sink_segment[nsink] = atoi(arg[8]); // # of segment for polygon sinks
  //   sink_mfp[nsink] = atof(arg[9]); // mean free path in this sink
  //   nabsorption[nsink] = 0; // initialize number of absorptions
  //   sink_creation(nsink); //create the nth sink, can overlap with other sinks
  //   nsink ++;
  // }
  else if (strcmp(command, "sink") ==0) {
    if(narg != 8) error->all(FLERR,"illegal sink command");
    if(sink_flag == 0) memory->create(isink, nlattice,"app/seg:isink");
    if(sink_flag == 0) {for(i = 0; i < nlattice; i++)  isink[i] = 0;} // set no sink initially
    sink_flag = 1;

    nsink ++;  // sink id starts from 1
    grow_sinks();
    sink_range[nsink] = atof(arg[0]); // thickness of sink
    sink_shape[nsink] = atoi(arg[1]); // 1 dislocation, 2 interface, 3 3D region
    xsink[nsink][0] = atof(arg[2]); // coordinaiton of sink center
    xsink[nsink][1] = atof(arg[3]);
    xsink[nsink][2] = atof(arg[4]);
    sink_normal[nsink] = atoi(arg[5]); // normal of planar sinks
    sink_radius[nsink] = atof(arg[6]); // radius for circular or polygonal sinks
    sink_segment[nsink] = atoi(arg[7]); // # of segment for polygon sinks
    sink_creation(nsink); //create the nth sink, can overlap with other sinks
  }

  // define element-sink interaction by eisink. One element and one sink in each line. The id of a sink is the order in sink commands, starting from 1
  else if (strcmp(command, "sink_interaction") ==0) {

    if(sink_flag == 0) error->all(FLERR,"sink_interaction needs to be defined after the sink command");
    if (narg < 3) error->all(FLERR,"Illegal sink_interaction command");
    if (eisink_flag == 0) {// create element-sink interaction.
       eisink_flag = 1;
       memory->create(eisink,nelement,nsink+1,"app/ris:eisink");
       memory->create(nabsorption,nelement,nsink+1,"app/ris:nabsorption");
       memory->create(nreserve,nelement,nsink+1,"app/ris:nreserve");
       memory->create(sink_mfp,nelement,nsink+1,"app/ris:sink_mfp");
       for (i = 0; i < nelement; i++ ) {
           for (j = 0; j < nsink+1; j++ ) {
               eisink[i][j] = 0.0;  // eisink = 0.0: no inteaction; eisink< -100: complete absortion; others: trapping or repulsion
               nabsorption[i][j] = 0;
               nreserve[i][j] = 0;
               sink_mfp[i][j] = 1.0;
           }
       }
    }


    int eletype = atoi(arg[0]); // element starts from 0
    int sinkid = atoi(arg[1]); // sink id starts from 1
    eisink[eletype][sinkid] = atof(arg[2]);
    if(narg > 3) sink_mfp[eletype][sinkid] = atof(arg[3]);
  }

  // define sinks to defects, could be dislocations or interfaces or 3D regions
  // else if (strcmp(command, "sink_motion") ==0) {
  //
  //   if(sink_flag == 0) error->all(FLERR,"sink_motion needs to be defined after the sink command");
  //   if(narg < 3) error->all(FLERR,"illegal sink_motion command");
  //   if(sinkmotion_flag == 0)  memory->create(sink_dr, nsink+1, "app/seg:sink_dr");
  //   if(sinkmotion_flag == 0)  memory->create(sink_dt, nsink+1, "app/seg:sink_dt");
  //   if(sinkmotion_flag == 0)  memory->create(sink_dt_new, nsink+1, "app/seg:sink_dt_new");
  //   if(sinkmotion_flag == 0)  memory->create(sink_dt_old, nsink+1, "app/seg:sink_dt_old");
  //   sinkmotion_flag = 1;
  //
  //   for (i = 0; i < nsink+1; i++) {sink_dr[i] = -1.0; sink_dt[i] = 0.0;} // by default no sink motion
  //
  //   int nseg = narg/3;
  //   for (i = 0; i < nseg; i++) {
  //       int sinkid = atoi(arg[i*3]); // sinkid start from 1
  //       if(sinkid > nsink) error->all(FLERR,"sink_motion needs to be defined after the sink command");
  //       sink_dr[sinkid] = atoi(arg[i*3+1]); // direction of sink motion (1 for x, 2 for y, and 3 or z), each step the displacement is a0/2
  //       double velocity = atof(arg[i*3+2]); // a0/s
  //       sink_dt[sinkid] = 0.5e12/velocity; // timestep for sink motion, in picosecond; the displacement each step is 0.5
  //   }
  // }
  // reactions for absorption and emission
  else if (strcmp(command, "reaction") ==0) {

    if(narg != 6) error->all(FLERR,"illegal reaction command");
    //if(diffusionflag) error->warning(FLERR,"MSD calculated with reactions"); //LC comment for current
    reaction_flag = 1;
    grow_reactions(); // grow reation list

    rsite[nreaction] = atoi(arg[0]); // reaciton site: type of lattice site
    rinput[nreaction] = atoi(arg[1]); // input element
    routput[nreaction] = atoi(arg[2]); // output element
    rbarrier[nreaction] = atof(arg[3]); // reaction barrier
    rrate[nreaction] = atof(arg[4]); // reaction rate
    rtarget[nreaction] = atoi(arg[5]); // target number of output element

    nreaction ++;
  }

  /* comment out, not using anymore, but keep it
  else if (strcmp(command, "surfacebarrier") ==0) {

    if (narg < 2 || narg % 2 != 0) error->all(FLERR,"Illegal surfacebarrier command");
    barrierflag = 1;
    memory->create(surfbarrier,nelement+1,"app/coros:surfbarrier");

    for (i=0; i<narg-1; i++) {
      if(i % 2 == 0){ j = atoi(arg[i]);
        surfbarrier[j] = atof(arg[i+1]);
      }
    }
  }
  */
  // command for initial salt potential, site and potential(double) not used now
  /*
  else if (strcmp(command, "saltpotential") ==0) {
    int saltsite;
    double saltpotential;
    if(narg != 2) error->all(FLERR,"illegal reaction command");
    saltsite = atoi(arg[0]); // salt site
    saltpotential = atoi(arg[1]); // initial potential
  }
  */
  // command for salt diffusion
  else if (strcmp(command, "saltdiffusion") ==0) {
    if(narg != 1) error->all(FLERR,"illegal salltdiffusion command");
    saltdiffusion_flag = 1;
    grow_saltdiffusion();
    salt_bfreq[nsaltdiffusion] = atoi(arg[0]); // diffusion frequency
    nsaltdiffusion ++; // number of saltdiffusion events
  }
  // command to enable propensity/barrier extraction
  else if (strcmp(command, "barrier_extract") ==0) {
    if(narg != 1) error->all(FLERR,"illegal barrier_extract command");
    extract_flag = atoi(arg[0]); // extract flag is on/enable if = 1, disable for other value
  }

  else if (strcmp(command, "evap_extract") ==0) {
    if(narg != 1) error->all(FLERR,"illegal evap_extract command");
    evap_extract_flag = atoi(arg[0]); // extract flag is on/enable if = 1, disable for other value
  }
  else if (strcmp(command, "np_extract") ==0) {
    if(narg != 1) error->all(FLERR,"illegal np_extract command");
    np_extract_flag = atoi(arg[0]); // np_extract flag is on/enable if = 1, disable for other value
  }
  else if (strcmp(command, "large_propensity_extract") ==0) {
    if(narg != 1) error->all(FLERR,"illegal large_propensity_extract command");
    large_propensity_extract_flag = atoi(arg[0]); // large_propensity_extract is on/enable if = 1, disable for other value
  }

  else if (strcmp(command, "dump_event") ==0) {
    if(narg != 1) error->all(FLERR,"illegal dump_event command");
    if (atoi(arg[0])==1){
      dump_event_flag = 1; //on/enable if = 1, disable for other value
    }
  }

  else if (strcmp(command, "dump_event_all") ==0) {
    if(narg != 1) error->all(FLERR,"illegal dump_event_all command");
    dump_event_all_flag = atoi(arg[0]); // dump_event_flag is on/enable if = 1, disable for other value
  }

  // stop function based on number of corrosion event
  else if (strcmp(command, "coros_stop") ==0) {

    if(narg != 1) error->all(FLERR,"illegal coros_stop command");
    threshold_Cr = atoi(arg[0]);
    KMC_stop_flag = 1;

  }

// take surface diffusion as exponential function
else if (strcmp(command, "surfaceeffect") ==0) {

  if(narg != 1) error->all(FLERR,"illegal surfaceeffect command");
  surface_effect_b = atoi(arg[0]);
}
    // feature for timming each functioncall
  else if (strcmp(command, "time_check") ==0) {
    if(narg != 1) error->all(FLERR,"illegal time_check command");
    if (atoi(arg[0])==1){
      time_check_flag = 1;
    }

  }
  // corrosion_info_output input setting
  else if (strcmp(command, "corrosion_info_output") ==0) {

    if(narg != 1) error->all(FLERR,"illegal surfaceeffect command");
    if(atoi(arg[0])==1){  // 1 = on, other = off
      coros_flag = 1;
      memory->create(coros_time,nlocal,"app/coros:coros_time");
      memory->create(coros_id,nlocal,"app/coros:coros_id");
    }
    //coros_flag = atoi(arg[0]);

  }
/*
  else if (strcmp(command, "ballistic") ==0) {
    if(narg < 2) error->all(FLERR,"illegal ballistic command");
    ballistic_flag = 1;
    grow_ballistic();
    bfreq[nballistic] = atoi(arg[0]); // mix rate
    rdamp[nballistic] = atof(arg[1]); // mix range
    if(min_bfreq > bfreq[nballistic]) min_bfreq = bfreq[nballistic];
    nballistic ++; // number of mixing events
  }
*/
  else error->all(FLERR,"Unrecognized command");
}

/* ----------------------------------------------------------------------
   define 2NN list based on 1NN list. In BCC lattice, each 2NN of site i is shared
   by 4 1NNs of site i as their 1NNs.
------------------------------------------------------------------------- */

void Appcoros::define_2NN()
{ int candidate[144],frequency[144];
  int i,j,k,jd,kd,n1nn,n2nn,njnn,ncandidate;

  memory->create(numneigh2,nlocal+nghost,"app, numneigh2");
  memory->create(neighbor2,nlocal+nghost,MAX2NN,"app, neighbor2");

  for (i = 0; i < nlocal+nghost; i++) {
    for (j = 0; j < 144; j++) candidate[j] = 0;
    for (j = 0; j < 144; j++) frequency[j] = 0;

    ncandidate = 0;
    n1nn = numneigh[i];
    n2nn = 0;
    for (j =0; j < n1nn; j++)
    {
      jd = neighbor[i][j];
      njnn = numneigh[jd];
      for (k = 0; k < njnn; k++)
      {
        kd = neighbor[jd][k];

        if (kd != i) {
          candidate[ncandidate] = kd;
          ncandidate++;
        }
      }
    }

    for (j = 0; j < ncandidate; j++) {
      if(frequency[j] < 0) continue; // jd already observed earlier
      jd = candidate[j];
      frequency[j]++;
      for (k = j+1; k < ncandidate; k++) {
        kd = candidate[k];
        if(kd == jd) {frequency[j]++;  frequency[k] = -1;}
      }
      if(frequency[j] == 4) {
        int ifbreak = 0;
  	for (int m = 0; m < n1nn; m ++) { // To screen 1NNs for fcc crystals
	    if(jd == neighbor[i][m]) ifbreak = 1;
        }
        if (ifbreak) continue;

        if (n2nn == MAX2NN) error->all(FLERR, "Two many 2nd nearest neighbors defined, please expand the simulation cell dimensions or MAX2NN");
        neighbor2[i][n2nn] = jd;
        n2nn++; // selected if shared by 4 NNs
      }
    }
    numneigh2[i] = n2nn;
  }
}

/* ----------------------------------------------------------------------
   set site value ptrs each time iarray/darray are reallocated
------------------------------------------------------------------------- */

void Appcoros::grow_app()
{
  type = iarray[0];   // lattice type; i1 in input
  element = iarray[1];  // element type; i2 in input
  potential = iarray[2]; // i3 contain potential of salt by LC d1 in input
  if(diffusionflag) {
    aid = iarray[3]; // initially set as global ID, must use set i4 unique in command line
    for (int i=0 ; i<nlocal;i++){
      aid[i] = i;
    }
    disp = darray; // msd; zero initially
  }
}

/* ----------------------------------------------------------------------
   define virtual function site_energy(int)
------------------------------------------------------------------------- */

double Appcoros::site_energy(int)
{
  return 0.0;
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */

void Appcoros::init_app()
{
  int i,j;

  // define second nearest neighbor
  if (engstyle ==2) define_2NN();

  if (firsttime) {
    firsttime = 0;

    echeck = new int[nlocal];
    memory->create(firstevent,nlocal,"app:firstevent");

    // esites must be large enough for 3 sites and their 1st neighbors
    esites = new int[3 + 3*maxneigh];
  }

  // site validity
  int flag = 0;
  for ( i = 0; i < nelement; i++) nsites_local[i] = 0;
//fprintf(screen, "nlocal:%d\n", nlocal); //LC test
  for ( i = 0; i < nlocal; i++) {
    //fprintf(screen, "type[i]:%d, element[i]:%d\n", type[i], element[i]); // LC test
    if (type[i] < BCC || type[i] > SALT) {flag = 1;}
    if (element[i] < FE || element[i] > SIA) flag = 1;
    nsites_local[element[i]-1]++;
  }

  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall) error->all(FLERR,"One or more sites have invalid values");

  // check if reactions need to be enabled or disabled

  if(reaction_flag) {
    for( i = 0; i< nreaction; i++) {
      rcount[i] = 0;
      renable[i] = 0;
      target_global[i] = 0;
      target_local[i] = 0;
    }
  }
  // initialize the time_list for salt diffusion
  if(saltdiffusion_flag) {
    for(i = 0; i < nsaltdiffusion; i ++) {
       salt_time_old[i] = 0;
       salt_time_new[i] = 0;
    }
  }
/*
  // initialize the time_list for ballistic mixing
  if(ballistic_flag) {
    for(i = 0; i < nballistic; i ++) {
       time_old[i] = 0;
       time_new[i] = 0;
    }
  }

  // initialize the time_list for ballistic mixing
  if(diffusionflag) {
    for(i = 0; i < ndiff; i ++) {
      for(j = 0; j < nlocal; j++){
         disp[i][j] = 0.0;
      }
    }
  }

  // initialize the time_list for ballistic mixing
  if(concentrationflag) {
    for(i = ndiffusion; i < ndiffusion + concentrationflag; i ++) {
      for(j = 0; j < nlocal; j++){
         disp[i][j] = 0.0;
      }
    }
  }
  // initiate parameters for vacancy trapping
  if(time_flag) {
    double nprcs = (double)(domain->nprocs);
    fvt = 1.0/nprcs;
    itrap = 0;
    itime_current = itime_old = 0;
    treal_me = takmc_me = 0.0;
    dt_real = dt_akmc = 0.0;
  }
*/

//initialize the concentration vectors
if(concentrationflag) {
 //concentration_field();
 dt_new = 0.0;
 dt_site_c_new = 0.0;
 for(i = 1; i <= nelement; i++) {
    ct[i] = 0.0;
    ct_new[i] = 0.0;

    //initiallize initial concentration
    for(j = 0; j < nlocal; j++) {
 ct_site[i][j] = 0.0;
 ct_site_new[i][j] = 0.0;
 // diag version
 ct_site_diag[i][j] = 0.0;
  ct_site_new_diag[i][j] = 0.0;
 //disp[i+ndiff][j] = 0.0;
    }
 }

//LC initiallize initial temp time and i2
time_sector = 0.0;
time_sector_diag =0.0;
  for(j = 0; j < nlocal; j++) {
      ct_site_temp_t[j] = 0.0;
      ct_site_temp_i2[j] = element[j];
      site_time[j] = 0.0;
      ct_site_temp_t_diag[j] = 0.0;
      ct_site_temp_i2_diag[j] = element[j];
      site_time_diag[j] = 0.0;

}
  for (j=0; j<8; j++){
site_time_interval[j] = 0.0;
}
 //initialize i3_site value
dt_i3_site_new =0.0;
 for(int i = 0; i < 5; i++) {
   for(int j = 0; j < nlocal; j++){
     i3_site[i][j] = 0.0;
     i3_site_new[i][j] = 0.0;
   }
 }

}

//initialize the concentration vectors for non concentration flag
// if(concentrationflag==0) {
//   for(i = 1; i <= nelement; i++) {
//      ct[i] = 0.0;
//      //initiallize initial concentration
//      for(j = 0; j < nlocal; j++) {
//   ct_site[i][j] = 0.0;
//      }
//   }
//   //initialize i3_site value
//   for(int i = 0; i < 5; i++) {
//     for(int j = 0; j < nlocal; j++){
//       i3_site[i][j] = 0.0;
//     }
//   }
//
// }


// cacculate sink statistics; LC
//sink_statistics();

if(coros_flag) {  // LC coros info output
  for( i = 0; i< nlocal; i++) {
    coros_id[i] = -1;
    coros_time[i] = 0.0;
  }
}


// bulk_vac initialization by LC
int z;
int temp_z=-1;
z_height=0;
xy_sum=0;
for (i = 0; i < nlocal; i++)
{
  // z position: xyz[i][2]
  z = int(xyz[i][2]*2);  // transfer z axis to 1D profile id

  if (z == 0){xy_sum ++;}  // sum of #site per cross-section
  if (z > temp_z) {temp_z = z;
    z_height ++;
  }
}
//fprintf(screen, "temp_z:%d, z_sum:%d, sum:%d\n",temp_z,z_height, sum );
memory->create(zprofile,z_height,"app:zprofile");
// 1D flatten: reset
for (i = 0; i < z_height; i++) // reset profile
{
  zprofile[i] = 0.0;
}

}

/* ---------------------------------------------------------------------- */
/* ---------------------------initiation of looing-up list -------------- */

void Appcoros::setup_app()
{
  for (int i = 0; i < nlocal; i++) echeck[i] = 0;

  // clear event list
  nevents = 0;
  for (int i = 0; i < nlocal; i++) firstevent[i] = -1;
  for (int i = 0; i < maxevent; i++) events[i].next = i+1;
  freeevent = 0;

  // set propensities from rates
  if(temperature == 0.0) error->all(FLERR,"Temperature cannot be 0.0 for app coros");
  if(nn1flag == 0) error->all(FLERR, "First neareast neighbor bond not defined: Appcoros");
  if(nn2flag == 0 && engstyle == 2) error->all(FLERR, "Second nearest neighbor bond not defined: Appcoros");
  if(barrierflag == 0) error->warning(FLERR, "Diffusion barrier not defined: Appcoros");

  double KB = 0.00008617;
  KBT = temperature * KB;
  // count number of each particle at initial structure
  int kd;
  for (int i = 0; i < nlocal; i++) {
    if(element[i] == FE) total_Ni++ ;
    if(element[i] == VACANCY) total_vac++ ;
    if(element[i] == CU) total_Cr++ ;

    //count surface atom
    int sum = 0;
    for (int n=0; n< numneigh[i]; n++){
      kd = neighbor[i][n];
      if (element[kd] == VACANCY) sum++;
    }
    if (element[i] == FE && sum>0){surface_Ni++;}
    if (element[i] == CU && sum>0){surface_Cr++;}
  }
  fprintf(screen,"surface_Ni %d surface_Cr %d\n",surface_Ni, surface_Cr );

  //count initial composition for metal by LC:
  double Ni_ratio = (double(total_Ni)/((double(total_Ni)+double(total_Cr))));
  double Cr_ratio = (double(total_Cr)/((double(total_Ni)+double(total_Cr))));

  // unit conversion, round up and convert to int
  Ni_ratio = int(Ni_ratio*100 + .5);
  Ni_ratio = (Ni_ratio)/100;
  Cr_ratio = int(Cr_ratio*100 + .5);
  Cr_ratio = (Cr_ratio)/100;

  Int_Ni = int(Ni_ratio *100);
  Int_Cr = int(Cr_ratio * 100);
  fprintf(screen,"initial Ni %d Cr %d \n",Int_Ni, Int_Cr );
}

/* ----------------------------------------------------------------------
   compute energy of site i
------------------------------------------------------------------------- */

double Appcoros::sites_energy(int i, int estyle)
{
  int j,jd,n1nn;
  double eng = 0.0;

  if(estyle == 0) return eng;
  //energy from 1NN bonds
  n1nn = numneigh[i];  //num of 1NN
  for (j = 0; j < n1nn; j++) {
    jd = neighbor[i][j];
    eng += ebond1[element[i]][element[jd]];
  }

  //energy from 2NN bonds
  if (estyle == 2) {
    int n2nn = numneigh2[i];
    for (j = 0; j < n2nn; j++) {
      jd = neighbor2[i][j];
      eng += ebond2[element[i]][element[jd]];
    }
  }
//!! changed by LC, the energy here is total bond energy for one site/atom,
//!! for some cases, the site energy is shared by neighbor and need to be divided by 2
  return eng;
}

/* ----------------------------------------------------------------------
   compute elastic energy
------------------------------------------------------------------------- */
/* elastic
double Appcoros::elastic_energy(int i, int itype)
{
  double eng = 0.0;
  double pressure = 0.0;
  pressure = -(stress[i][0] + stress[i][1] + stress[i][2])/3.0;
  eng = pressure * evol[itype];
  return eng;
}
*/
/* ----------------------------------------------------------------------
  compute barriers for an exchange event between i & j
------------------------------------------------------------------------- */

double Appcoros::site_SP_energy(int i, int j, int estyle, double s)
{
  double eng = 0.0;
  double eng0i, eng0j, eng1i, eng1j; //energy before and after jump
  int iele = element[i];
  int jele = element[j];
  eng0i = sites_energy(i,estyle); //broken bond with i initially,
  eng0j = sites_energy(j,estyle); //broken bond with j initially
  // switch the element and recalculate the site energy
  element[i] = jele;
  element[j] = iele;
  eng1i = sites_energy(i,estyle); //broken bond with i initially,
  eng1j = sites_energy(j,estyle); //broken bond with j initially
  // switch back
  element[j] = jele;
  element[i] = iele;
  //barrier = migbarrier * surface_effect + (eng_after - eng_before)/2.0;
  eng = mbarrier[element[j]] * s + (eng1i + eng1j - eng0i -eng0j)/2; //LC add surface effect

  //add elastic contribution if applicable
  /* comment because assume no elastic interaction
  if(elastic_flag) {
    double eng_ei = elastic_energy(j,iele) - elastic_energy(i,iele);
    double eng_ej = elastic_energy(i,jele) - elastic_energy(j,jele);
    eng += (eng_ei + eng_ej)/2.0;
  }
*/
  return eng;
}

/* ----------------------------------------------------------------------
   KMC method
   compute total propensity of owned site summed over possible events
------------------------------------------------------------------------- */

double Appcoros::site_propensity(int i)
{
  int j, iid, jid;
  int k, kd, sum_kd;

  // !!parameter for surface diffusion by LC
  int sum_bond,sum_total, r;
  int num_iter = 0;
  int sum_atom = 0;
  double bond_ratio;
  double b;
  double surface_effect;
  // valid hop and recombination tabulated lists
  // propensity for each event is input by user
  clear_events(i);
  double prob_reaction = 0.0;
  double prob_hop = 0.0;
  double ebarrier = 0.0;
  double hpropensity = 0.0;

  //LC test
  // this test print the cartesion position of site i
  //fprintf(screen, "site[i]:%d, element[i]:%d\n", i, element[i]);
  //double **xyz = app->xyz;
  //double x = app->xyz[0][0];
  //double y = app->xyz[0][1];
  //double z = app->xyz[0][2];
  //xyz[clist[i]][0]
  //fprintf(screen, "x,y,z: %f, %f, %f\n", xyz[i][0],xyz[i][1],xyz[i][2]);
  //LC test



  // propensity for reactions, only when flagged and enabled
  // zero barrier event are dealted with separately in site_events()
  // barriers are from input
  if(reaction_flag) { //reaction flag
    for(j = 0; j < nreaction; j++) {
      if(renable[j] == 0) continue;
      if(element[i] == rinput[j] && type[i] == rsite[j]) {
        iid = rinput[j];
        jid = routput[j];

        // count number of neighbor salt by LC
        sum_kd = 0;
        for (k = 0; k < numneigh[i] ; k++){
          kd = neighbor[i][k];
          if (potential[kd] == 1){
            sum_kd ++;
          }
        }
        //if (sum_kd == 0) return prob_reaction;
        if (sum_kd == 0) break; // no reaction // edited may need to double check later
        // count number of neighbor salt

        // LC comment
        //if(!(sink_flag && isink[i][jid -1] == 1)) {//production at sinks allowed
          ebarrier = rbarrier[j];
          // elastic
          /*
          if(elastic_flag) ebarrier += elastic_energy(i,jid) - elastic_energy(i,iid);
          */
          // k = rate * concentration * exp(-Ea/KbT)
          //LC test comment
          hpropensity = rrate[j] * sum_kd * exp(-ebarrier/KBT);

          //LC test
          //hpropensity = rrate[j];

          // barrier/propensity print out by LC !!need double check
          // if (extract_flag==1){
          //   surface_effect = 1;
          //   r = 2;
          //   barrier_print(r, hpropensity, rrate[j], surface_effect, ebarrier);
          // }
          //
          //fprintf(screen,"reaction_propensity: %e\n",hpropensity);

          // extract large propensity case info for propensity > 1
          r =2; //rstyle ==2 for reaction
          if (large_propensity_extract_flag ==1 && hpropensity >= 1.0){
            //np_check(i,jid);
            barrier_print(r, hpropensity, attemptfrequency[element[jid]], surface_effect, ebarrier);
          }
          // extract large propensity case info for propensity > 1

          add_event(i,jid,2,j,hpropensity);
          prob_reaction += hpropensity;
        //}
      }
    }
  }

  if (element[i] != VACANCY) return prob_reaction;
  // check if acceleration is needed and update propensity if so
  /*
  if (acceleration_flag == 1 ) {
     for (j = 0; j < numneigh[i]; j++) {
       jid = neighbor[i][j];
       for ( int k = 0; k < ntrap; k ++) {
          if (element[jid] == trap_type[k]) return add_acceleration_event(i,jid);
       }
     }
  }
*/
  // for hop events, vacancy only currently
  // propensity calculated in site_SP_energy();
  for (j = 0; j < numneigh[i]; j++) {
    jid = neighbor[i][j];
    if(element[jid] != VACANCY) { // no vacancy-vacancy switch
      sum_atom++;
      // surface effect section
      sum_bond = 0;
      sum_total = 0;
      bond_ratio = 0.0;
      for (k = 0; k < numneigh[jid]; k++){
        sum_total++;
        if (element[neighbor[jid][k]]!= VACANCY){
            sum_bond++;
        }
      }
      bond_ratio = (sum_bond+1.0)/sum_total;
      r = 1; // rstyle = 1
      // surface effect = exp(-(NTotal- (Nbond +1 ))/b)
      surface_effect = 0.0;
      surface_effect = exp(-(sum_total-(sum_bond+1.0))/surface_effect_b);
      //
      ebarrier = site_SP_energy(i,jid,engstyle, surface_effect); // diffusion barrier
      //LC comment
      hpropensity = attemptfrequency[element[jid]] * exp(-ebarrier/KBT);

      //fprintf(screen,"hpropensity:%e\n", hpropensity);
      //LC test:  1.43*10**-7
      //hpropensity = 0.000000143;

      // extract all case info
      if (extract_flag==1){
          np_check(i,jid);
          barrier_print(r, hpropensity, attemptfrequency[element[jid]], surface_effect, ebarrier);
      }

      // extract non-positive barrier case info
      if (np_extract_flag ==1 && ebarrier <= 0&& r==1){
        np_check(i,jid);
        barrier_print(r, hpropensity, attemptfrequency[element[jid]], surface_effect, ebarrier);
      }

      // extract large propensity case info for propensity > 1
      if (large_propensity_extract_flag == 1 && hpropensity >= 1.0){
        barrier_print(r, hpropensity, attemptfrequency[element[jid]], surface_effect, ebarrier);
      }
      // extract large propensity case info for propensity > 1

      add_event(i,jid,1,-1,hpropensity);
      prob_hop += hpropensity;

      //LC check
        // if (element[jid] == element[i]){
        //   fprintf(screen,"i:%d, element[i]:%d, jid:%d, element[j]:%d propensity:%e \n",i, element[i], jid, element[jid], hpropensity);
        // }

       // if (jid==9){
       //   fprintf(screen, "check point:i:%d, element[jid]:%d\n", i, element[jid]);
       // }
    }
  }
  // extract evaporation case info
  if (evap_extract_flag ==1 && sum_atom == 1){
    np_check(i,jid);
    barrier_print(r, hpropensity, attemptfrequency[element[jid]], surface_effect, ebarrier);
  }

  return prob_hop + prob_reaction;
}

/* ----------------------------------------------------------------------
   KMC method
   choose and perform an event for site
------------------------------------------------------------------------- */

void Appcoros::site_event(int i, class RandomPark *random)
{
  int j,k,l,m,n,ii;

  // perform events with non_zero barriers
  // pick one event from total propensity by accumulating its probability
  // compare prob to threshhold, break when reach it to select event
  double threshhold = random->uniform() * propensity[i2site[i]];
  double proball = 0.0;

  // check if trapped by a solute
  if(time_flag) itrap = vacancy_trap(i);
  // find the event to perform
  int ievent = firstevent[i];
  while (1) {// LC note: select event based on event group of site i
    proball += events[ievent].propensity;
    //if(element[events[ievent].jpartner]==2){fprintf(screen,"check \n");} // LC check
    if (proball >= threshhold) break;
    ievent = events[ievent].next;
  }

  // perform hop or reaction event
  int rstyle = events[ievent].style;
  int which = events[ievent].which; // type of reactions or neighbor id for acceleration
  j = events[ievent].jpartner;

//LC check
  // if (element[j] == element[i]){
  //   fprintf(screen,"ievent: %d,i:%d, element[i]:%d, jid:%d, element[j]:%d propensity:%e \n",ievent,i, element[i], j, element[j], events[ievent].propensity);
  // }


  int kp, lp;   //LC swap salt potential also, if swap with vac with salt
  // switch element between site i and jpartner for hop diffusion
  if(rstyle == 1 || rstyle == 3 || rstyle ==4) {
    //k = element[i];
    //kp = potential[i]; // LC swap salt potential
    if (rstyle == 4) { // switch with a 1NN of which
       k = element[i];
       kp = potential[i]; // LC swap salt potential
       l = element[which];
       lp = potential[which]; //LC swap salt potential
       element[i] = l;
       potential[i] = lp; //LC swap salt potential
       element[which] = element[j];
       potential[which] = potential[j]; //LC swap salt potential
       element[j] = k;
       potential[j] = kp; // LC swap salt potential
    } else { // switch with a 1NN of i
      k = element[i];
      kp = potential[i]; // LC swap salt potential
      element[i] = element[j];
      element[j] = k;
      potential[i] = potential[j] ; //LC swap salt potential
      potential[j] = kp; //LC swap salt potential
    }

//LC check
    if (element[j] == element[i]){
      fprintf(screen,"ievent: %d,i:%d, element[i]:%d, jid:%d, element[j]:%d propensity:%e \n",ievent,i, element[i], j, element[j], events[ievent].propensity);
    }

    update_propensity(i);
    update_propensity(j);

    //LC update current i and j
    cur_i = i;
    cur_j = j;
    // int jid;
    // for (int idx=0;idx<numneigh[i];idx++){
    //   jid = neighbor[i][idx];
    //   update_propensity(jid);
    // }
    // for (int idx=0;idx<numneigh[j];idx++){
    //   jid = neighbor[j][idx];
    //   update_propensity(jid);
    // }
    //update 2NN of i and j

    hcount[element[i]] ++; // hop count

    // this part is count number of diffusion for each elements by LC
    // not use but keep it
    /*
    if(element[i] == 1 && element[j] == 2){ // bulk diff of id2 = 1
      nbulkfe ++;
    }
    if(element[i] == 3 && element[j] == 2){ // bulk diff of id2 = 3
      nbulkcu ++;
    }
    */

    // calculate MSD for each atom if activated
    if(diffusionflag) {
      // switch global atomic id
      k = aid[i];
      aid[i] = aid[j];
      aid[j] = aid[i];

      //update and switch displacement
      int periodicity[3];
      double dij[3],lprd[3];
      periodicity[0] = domain->xperiodic;
      periodicity[1] = domain->yperiodic;
      periodicity[2] = domain->zperiodic;
      lprd[0] = domain->xprd;
      lprd[1] = domain->yprd;
      lprd[2] = domain->zprd;
//LC note: if go across boundary as positive, dij<0 and dij~lprd[k] (length of axis length)
      for (k = 0; k < 3; k++) { //update
        dij[k] = xyz[j][k] - xyz[i][k];
        if (periodicity[k] && dij[k] >= lprd[k]/2.0) dij[k] -= lprd[k];
        if (periodicity[k] && dij[k] <= -lprd[k]/2.0) dij[k] += lprd[k];
        disp[k][i] += dij[k];
        disp[k][j] -= dij[k];
      }

      for (k = 0; k < 3; k++) { //switch
         dij[k] = disp[k][i];
         disp[k][i] = disp[k][j];
         disp[k][j] = dij[k];
      }
      disp[3][i] = disp[0][i]*disp[0][i] + disp[1][i]*disp[1][i] + disp[2][i]*disp[2][i];
      disp[3][j] = disp[0][j]*disp[0][j] + disp[1][j]*disp[1][j] + disp[2][j]*disp[2][j];


    }

  } else {// LC note: rstyle ==2 in this else session, and for reaction only

 // reaction events with non-zero barriers
    k = element[i];
    element[i] = j;
    rcount[which] ++;
    nsites_local[k-1] --;
    nsites_local[j-1] ++;

    salt_remove(i);


    if(coros_flag) corrosion_info_output(i);  // LC Cr info record, output at the end

    total_vac ++;  // real time number of i2
    total_Cr --;
    //update after reaction, it will update in salt_remove function
    //just leave here for check
    update_propensity(i);
//update_propensity(j);
    // update reaction target number
    for(ii = 0; ii < nreaction; ii++) {
      if(routput[ii] != j) continue;
      target_local[ii] --;
    }
  }



  // perform zero_barrier events: absorption and recombination
  if(rstyle == 1 || rstyle == 3 || rstyle == 4) {

    //sink absorption,for hop only since no element produced at its sinks
    if(nsink > 0) {
      absorption(i);
      absorption(j);
    }
    update_propensity(i);
    update_propensity(j);
// LC comment
    // if(nsink > 0) {
    //   k = element[j];
    //   for (n = 0; n < nsink; n++) {
    //     if(isink[j][k-1] == 1) {
    //       double rand_me = rancoros->uniform();
    //       if(rand_me <= 1.0/sink_mfp[n]) {
    //         nabsorption[n] ++;
    //         nsites_local[k-1] --;
    //         element[j] = FE;
    //         // update reaction target number
    //         if(reaction_flag == 1) {
    //           for(ii = 0; ii < nreaction; ii++) {
    //             if(routput[ii] != k) continue;
    //             target_local[ii] ++;
    //           }
    //         }
    //       }
    //     }
    //   }
    // }

  //check recombinations after hopping; recombine with first SIA found
  /*
    for (n = 0; n < numneigh[j]; n++) {
      m = neighbor[j][n];
      if(element[j] == VACANCY && element[m] == SIA) {
        nrecombine ++;
        nsites_local[VACANCY-1] --;
        nsites_local[SIA-1] --;
        element[j] = FE;
        element[m] = FE;
        // update reaction target number
        if(reaction_flag == 1) {
          for(ii = 0; ii < nreaction; ii++) {
            if(routput[ii] != VACANCY && routput[ii] != SIA) continue;
            target_local[ii] ++;
          }
        }
      }
    }
*/
  } else {

    //for a vavancy produced (rstyle 2), check if recombine with neighbor SIAs
    /*
    for (n = 0; n < numneigh[i]; n++) {
      m = neighbor[i][n];
      if(element[i] == VACANCY && element[m] == SIA) {
        nrecombine ++;
        nsites_local[VACANCY-1] --;
        nsites_local[SIA-1] --;
        element[i] = FE;
        element[m] = FE;
        // update reaction target number
        if(reaction_flag == 1) {
          for(ii = 0; ii < nreaction; ii++) {
            if(routput[ii] != VACANCY && routput[ii] != SIA) continue;
            target_local[ii] ++;
          }
        }
      }
    }
    */
  }


//!! disable interface check by LC
/*
  // update type and reassign metal, interface and salt regions
  // for diffusion check and reaction check
  // update_region include subclass update_region_loop
  update_region(i, j, rstyle);
  // compute propensity changes for participating sites i & j and their neighbors
  update_propensity(i);
  update_propensity(j);
*/

  // check if any active reactions needs to be disabled
  /*
  if(reaction_flag == 1) {
    for(ii = 0; ii < nreaction; ii ++) {
      if(target_local[ii] <= 0 && renable[ii] == 1) {
        renable[ii] = 0;
        reset_propensity();
        return;
      }
    }
  }
  */
}

/* ----------------------------------------------------------------------
   update propensity for site i,j and their neighbors after exchange
   ignore update of sites with i2site < 0; updated during each sweep
   use echeck[] to avoid resetting propensity of same site
------------------------------------------------------------------------- */
void Appcoros::update_propensity(int i)
{
  int m,n;
  int nsites = 0;
  int isite = i2site[i];

  //if (isite < 0) return;
  // LC edite
  // this update should update i, j and 1NN of i and j,
  // the previous code will not update 1NN of i if i is out of sector
  if(isite>=0){
    propensity[isite] = site_propensity(i);
    esites[nsites++] = isite;
    echeck[isite] = 1;
   }

  for (n = 0; n < numneigh[i]; n++) { // update propensity for 1NN of site i
    m = neighbor[i][n];
    isite = i2site[m];
    if (isite >= 0 && echeck[isite] == 0) {
      propensity[isite] = site_propensity(m);
      esites[nsites++] = isite;
      echeck[isite] = 1;
    }
  }

  if(engstyle == 2 ) {
    for (n = 0; n < numneigh2[i]; n++) { // update propensity for 2NN of site i
      m = neighbor2[i][n];
      isite = i2site[m];
      if(isite >= 0 && echeck[isite] == 0) {
        propensity[isite] = site_propensity(m);
        esites[nsites++] = isite;
        echeck[isite] = 1;
      }
    }
  }

  solve->update(nsites,esites,propensity);

  // clear echeck array
  for (m = 0; m < nsites; m++)
    echeck[esites[m]] = 0;
}

/* ----------------------------------------------------------------------
   check if any reaction is need to be enabled or disabled
------------------------------------------------------------------------- */

void Appcoros::check_reaction()
{
  int i,m,n,n_me,n_total,flag_local,nsites_global;
  double nprcs = (double)(domain->nprocs);

  //update the statistics of each elements locally
  for (i = 0; i < nelement; i++) nsites_local[i] = 0;
  for (i = 0; i < nlocal; i++) nsites_local[element[i]-1]++;

  //compute how many atoms need to be generated globally
  for (i = 0; i < nelement; i++) {
    nsites_global = 0;
    MPI_Allreduce(&nsites_local[i],&nsites_global,1,MPI_INT,MPI_SUM,world);

    for(n = 0; n < nreaction; n++){
      if(i+1 == routput[n]) target_global[n] = rtarget[n] - nsites_global;
    }
  }

  //assign the expected atoms evenly to local processors, turn on reaction locally if needed
  flag_local = 0;
  for(n = 0; n < nreaction; n++){
    target_local[n] = 0;
    m = renable[n];
    renable[n] = 0;

    if(target_global[n] > 0) {
      for( i = 0; i < target_global[n]; i++){
        n_total = 0;

        while(n_total < 1) {
          n_me = 0;
          double rand_me = rancoros->uniform();
          if(rand_me < 1.0/nprcs) n_me += 1;
          target_local[n] += n_me;
          MPI_Allreduce(&n_me,&n_total,1,MPI_INT,MPI_SUM,world);

          if(n_total > 1) { //only increase by one allowed in each round
            n_total = 0;
            target_local[n] -= n_me;
          }
        }
      }
    }

    if(target_local[n] > 0) renable[n] = 1;
    flag_local += abs(renable[n] - m);
  }

  // reset site propensity if any enabling or disabling occurs
  if(flag_local > 0) reset_propensity();


}

/* ----------------------------------------------------------------------
   if any reaction is enabled or disabled, reset the propensity
   for all local lattice sites
------------------------------------------------------------------------- */
void Appcoros::reset_propensity()
{
  for (int i = 0; i < nset; i++) { //need update all sectors
    for (int m = 0; m < set[i].nlocal; m++){
      set[i].propensity[m] = site_propensity(set[i].site2i[m]);
      set[i].solve->update(m,set[i].propensity);
    }
  }
}

/* ----------------------------------------------------------------------
   clear all events out of list for site I
   add cleared events to free list
------------------------------------------------------------------------- */

void Appcoros::clear_events(int i)
{
  int next;
  int index = firstevent[i];

  while (index >= 0) {
    next = events[index].next;
    events[index].next = freeevent;
    freeevent = index;
    nevents--;
    index = next;
  }
  firstevent[i] = -1;
}

/* ----------------------------------------------------------------------
  add an acceleration event to list for site I
  return a propensity
------------------------------------------------------------------------- */
/*
double Appcoros::add_acceleration_event(int i, int j)
{
  int ni,nj,nn,nid,njd;
  double pr,pl,texit,p1,p2,p3,p4,w1,w2,w3,w4,eb;
  ni = numneigh[i];
  nj = numneigh[j];
  double *hpi = new double[ni];
  double *hpj = new double[nj];
  w1 = w2 = w3 = w4 = 0.0;
  p1 = p2 = p3 = p4 = 0.0;
  nn = 0;
  w1 = exp(-site_SP_energy(i,j,engstyle)/KBT);
  for (int m = 0; m < ni; m ++) {
      nid = neighbor[i][m];
      if (nid == j) continue;
      if (element[nid] == VACANCY) { hpi[nn] = 0.0;
      } else {
        eb = site_SP_energy(i,nid,engstyle); // diffusion barrier
        hpi[nn] = exp(-eb/KBT);
      }
      w2 += hpi[nn];
      nn ++;
  }
  int Ei = element[i];
  int Ej = element[j];
  element[i] = Ej;
  element[j] = Ei;
  nn = 0;
  w3 = exp(-site_SP_energy(j,i,engstyle)/KBT);
  for (int n = 0; n < nj; n ++) {
      njd = neighbor[j][n];
      if (njd == i) continue;
      if (element[njd] == VACANCY) { hpj[nn] = 0.0;
      } else {
        eb = site_SP_energy(j,njd,engstyle); // diffusion barrier
        hpj[nn] = exp(-eb/KBT);
      }
      w4 += hpj[nn];
      nn ++;
  }
  p1 = w1 / (w1+w2);
  p2 = w2 / (w1+w2);
  p3 = w3 / (w3+w4);
  p4 = w4 / (w3+w4);
  double t1 = 1.0/(w1+w2);
  double t2 = 1.0/(w3+w4);
  texit = (p2*t1 + p2*t2*p1*p3 + p1*p4*t1 + p1*p4*t2)/(1-p1*p3)/(1-p1*p3);
  pl = p1 / (1-p1*p3);
  pr = p1*p4 / (1-p1*p3);
  // add site events
  nn = 0;
  for (int m = 0; m < ni; m ++) {
      nid = neighbor[i][m];
      if (nid == j) continue;
      hpi[nn] *= pl;
      hpi[nn] /= w2; // normalize the propensity
      add_event(i,nid,3,j,hpi[nn]); // rstyle == 3 means exit left
      nn ++;
  }
  nn = 0;
  for (int n = 0; n < nj; n ++) {
      njd = neighbor[j][n];
      if (njd == i) continue;
      hpj[nn] *= pr;
      hpj[nn] /= w4; // normalized the propensity
      add_event(i,njd,4,j,hpj[nn]); // rstyle == 4 means exit right
      nn ++;
  }
  element[i] = Ei;
  element[j] = Ej;
  return 1.0 / texit; // return the propensity
}
*/
/* ----------------------------------------------------------------------
  add an event to list for site I
  event = exchange with site J with probability = propensity
------------------------------------------------------------------------- */

void Appcoros::add_event(int i, int j, int rstyle, int which, double propensity)
{
  // grow event list and setup free list
  if (nevents == maxevent) {
    maxevent += DELTAEVENT;
    events =
      (Event *) memory->srealloc(events,maxevent*sizeof(Event),"app:events");
    for (int m = nevents; m < maxevent; m++) events[m].next = m+1;
    freeevent = nevents;
  }

  int next = events[freeevent].next;

  //for hop, rstyle = 1, and switch element at sites I and which
  //for reaction, rstyle = 2, and switch element at site I to which
  //events[freeevent].kpartner = kpartner;
  events[freeevent].style = rstyle;
  events[freeevent].jpartner = j;
  events[freeevent].which = which;
  events[freeevent].propensity = propensity;

  events[freeevent].next = firstevent[i];
  firstevent[i] = freeevent;
  freeevent = next;
  nevents++;

  // dump all eventlist to one file; !!need double check
  if(dump_event_all_flag ==1){
  dump_event_all(rstyle, which, j, next, propensity);}

}

/* ----------------------------------------------------------------------
   check if perform ballistic mixing
------------------------------------------------------------------------- */
void Appcoros::check_ballistic(double t)
{
  int nmix = 0;
  for(int i = 0; i < nballistic; i ++) {
     time_new[i] = static_cast<int>(t/bfreq[i]);
     nmix = time_new[i] - time_old[i];

     while (nmix) {  //perform mixing nmix times
       nmix --;
       ballistic(i);
       if(nmix == 0) time_old[i] = time_new[i];  //update time
    }
  }
}

/* ----------------------------------------------------------------------
  perform ballistic mixing
------------------------------------------------------------------------- */

void Appcoros::ballistic(int n)
{
  int i,iwhich,iid,jid,ibtype,jbtype,ibtype_all,jbtype_all;
  int iproc,jproc,found_me,found_all;
  double rand_me,xiid[3];

  // find atom j for the exchange
  iproc = jproc = -1;
  ibtype = jbtype = ibtype_all = jbtype_all = 0;
  xmix[n][0] = xmix[n][1] = xmix[n][2] = 0.0;
  xiid[0] = xiid[1] = xiid[2] = 0.0;

  found_me = found_all = 0;
  while (found_all == 0) {
    rand_me = rancoros->uniform();
    if(rand_me < static_cast<double> (1.0/nprocs)) found_me ++;
    MPI_Allreduce(&found_me, &found_all,1,MPI_INT,MPI_SUM,world);

    // only one processor can be selected each time
    if(found_all > 1) {
      found_me = 0;
      found_all = 0;
    }
  }

  iid = -1;
  if(found_me) {
    iproc = me;
    iwhich = -1;
    while(iwhich < 0 || iwhich >= nlocal) {
      rand_me = rancoros->uniform();
      iwhich = static_cast<int> (nlocal*rand_me);}

    iid = iwhich;
    ibtype = element[iid];
    xiid[0] = xyz[iid][0];
    xiid[1] = xyz[iid][1];
    xiid[2] = xyz[iid][2];
  }

  MPI_Allreduce(&xiid[0],&xmix[n][0],1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&xiid[1],&xmix[n][1],1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&xiid[2],&xmix[n][2],1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&ibtype,&ibtype_all,1,MPI_INT,MPI_SUM,world);

  // calculate the total exchange propensity
  ballistic_probability(n);

  // find proc j and then jid based on the location of iid
  found_me = found_all = 0;
  while (found_all == 0) {
    rand_me = rancoros->uniform();
    if(rand_me < pn_local[n]/pn_global[n]) found_me ++;
    MPI_Allreduce(&found_me, &found_all,1,MPI_INT,MPI_SUM,world);

    // only one processor can be selected each time
    if(found_all > 1) {
      found_me = 0;
      found_all = 0;
    }
  }

  // find atom j on selected processor
  jid = -1;
  if(found_me) {
    i = 0;

    rand_me = rancoros->uniform()*pn_local[n];
    while(jid < 0 && i < nlocal) {
      if(rand_me <= pmix[n][i]) jid = i;
      i ++;
    }

    if(jid < 0) error->all(FLERR,"could not find atom j for mixing !");

    jproc = me; // proc owns jid
    jbtype = element[jid];
  }

  MPI_Allreduce(&jbtype,&jbtype_all,1,MPI_INT,MPI_SUM,world);


  // exchange type of iid and jid and update propensity
  if(me == iproc) {
    element[iid] = jbtype_all;
    update_propensity(iid);
  }

  if(me == jproc) {
    element[jid] = ibtype_all;
    update_propensity(jid);
  }

}

/* ----------------------------------------------------------------------
   check if any sink motion needs to be performed
   // LC not use now
------------------------------------------------------------------------- */
void Appcoros::check_sinkmotion(double t)
{
  int nmove = 0;
  for(int i = 1; i < nsink+1; i ++) {
     if(sink_dr[i] < 0.0) continue; // skip static sinks
     sink_dt_new[i] = static_cast<int>(t/sink_dt[i]);
     nmove = sink_dt_new[i] - sink_dt_old[i];
     while (nmove > 0) {  //perform mixing nmix times
       nmove --;
       sink_motion(i);
       if(nmove == 0) sink_dt_old[i] = sink_dt_new[i];  //update time
    }
  }

  return;
}

/* ----------------------------------------------------------------------
   move sink n in the direction of sink_dr[n] with an amount of a0/2
   This is done by delete sink n and recreate it with the center displaced by a0/2
   // LC not use now
------------------------------------------------------------------------- */
void Appcoros::sink_motion(int n)
{
  int nlattice = nlocal + nghost;
  double dr=sink_dr[n];
  if(dr == 0.0) xsink[n][0] = xsink[n][0] + 0.5;
  if(dr == 1.0) xsink[n][1] = xsink[n][1] + 0.5;
  if(dr == 2.0) xsink[n][2] = xsink[n][2] + 0.5;

  for (int i=0; i<nlattice; i++) {
      if(isink[i] == n) isink[i] = 0; //remove sink n
  }

  sink_creation(n); //recreate sink n
  return;
}

/* ----------------------------------------------------------------------
   absorption of an element at site i. This is probabbly not needed.
   LC modify:
   the sink absorption remove vac and refill by Ni/Cr with chance 80% for Ni
   20% for Cr.
   previous code just used the ratio of current model, which has bias on possiblity
   I also commented out most of unused function which improve some of efficiency
------------------------------------------------------------------------- */
void Appcoros::absorption(int i)
{
  int j,k,m,n,ei,ejd,ntotal,rand_me,jbound,jid;
  int temp_element;

// for(int id = 1; id< 100; i++){
//   rand_me = static_cast<int>(rancoros->uniform()*100);
//   fprintf(screen, "rand_me: %d\n", rand_me);
// }
  n = isink[i]; //sink id
  if(n <= 0) return; // site i is not a sink

  ntotal = 0;
  ei = element[i];
  ejd = -1;

  if (ei!= VACANCY) return; // only absorb VAC
  //if(ei > INT) return; // currently only absorbs VAC and INT, LC comment
  if(eisink[ei][n] > -100) return; // not an absorbing sink
  //if(rancoros->uniform() > 1.0/sink_mfp[ei][n]) return; // no absorption occurs , LC comment

  nabsorption[ei][n] ++; // number of element ei absorbed by type of isink[i] sink
  if(ei == VACANCY)  { // for vacancy choose a reserved atom to occupy
    // LC comment
     //for(m = CE1; m < nelement ; m++) {if(nreserve[m][n] > 0) ntotal += nreserve[m][n];} // count all reserved SIAs at sinks
     // LC modify
     while(ejd < 0) {
         //rand_me = static_cast<int>(rancoros->uniform()*nlocal); // LC note: incorrect
         rand_me = static_cast<int>(rancoros->uniform()*100);
         //fprintf(screen, "rand_me: %d\n", rand_me);
         // change the rand_me setting for initial Ni composition and threshold
         if (rand_me <= Int_Ni){temp_element = 1;}
         else{temp_element = 3;}
         ejd = temp_element;
         //fprintf(screen, "ejd: %d\n", ejd);
         //if(element[rand_me] != VACANCY) ejd = element[rand_me];
     }
     // LC comment
     // if(ntotal > 0) {// choose from reserved elements at sinks
     //   j = CE1;
     //   jbound = 0;
     //   rand_me = static_cast<int>(rancoros->uniform()*ntotal);
     //   while(ejd < 0) {
	   //       if(nreserve[j][n] > 0) jbound += nreserve[j][n];
     //       if(rand_me < jbound) ejd = j;
	   //       j++;
	   //       if(j == nelement) ejd = j-1;
     //   }
     // } else {// if none reserved "borrow" one proportionally to the nominal composition from reserved
     //   while(ejd < 0) {
     //       rand_me = static_cast<int>(rancoros->uniform()*nlocal);
	   //       if(element[rand_me] != VACANCY) ejd = element[rand_me];
     //   }
     // }
     if(ejd == -1) error->all(FLERR,"no element found to recombine with vacancy at sink");

     // recombine and recount the # of atoms and reserved atomsat sink n
     nsites_local[ei] --;
     nsites_local[ejd] ++;
     element[i] = ejd;
     // check potential value for absorption
     if (element[i] == 1){potential[i] = 3;}
     if (element[i] == 3){potential[i] = 4;}
     //nreserve[ejd][n] --; // LC comment for segmentation fault
  } else if (ei != VACANCY) { // randomly reserve one element from a dumbbell
    // LC comment all, due to no dumbell in this code, and ei must be vac at this case
     // nsites_local[ei] --;
     //
     // // mix the two dumbbell atoms with the atoms in sink n
     // // then ramdomly select one to be put into reservior
     // // the probability of dumb1 or dumb2 or a sink atom to go to resevior is 1/(nsink_site+1)
     //
     // j = static_cast<int>(rancoros->uniform()*(sinksite[n-1]+2));
     // if(j == sinksite[n-1]+1) { //dmb2 goes to reservior
     //    element[i] = dmb1[i];
     //    nreserve[dmb2[i]][n] ++;
     // } else if(j == sinksite[n-1]) { //dmb1 goes to reservior
     //    element[i] = dmb2[i];
     //    nreserve[dmb1[i]][n] ++;
     // } else { // a selected sink atom goes to reservior, dmb1 stays at i and dmb2 goes to the sink site
     //    element[i] = dmb1[i];
	   //    jid = sinkid[n-1][j];
     //    nsites_local[element[jid]] --;
     //    nreserve[element[jid]][n] ++;
     //
     //    element[jid] = dmb2[i];
     //    nsites_local[dmb2[i]] ++;
     // }
     //
     // nsites_local[element[i]] ++;
     // dmb1[i] = -1;
     // dmb2[i] = -1;
     // siatype[i] = -1;
  }

  // update reaction target number
  // LC comment out for not using this section
  /*
  if (reaction_flag == 1) {
     for(j = 0; j < nreaction; j++) {
        if(routput[j] == ei) target_local[j] ++; // one less ei atom
        if(routput[j] == element[i]) target_local[j] --; // one more element[i] atom

        if(target_local[j] <= 0 && renable[j] == 1) {
          renable[j] = 0;
          reset_propensity();
        }
        if(target_local[j] > 0 && renable[j] == 0) {
          renable[j] = 1;
          reset_propensity();
        }
     }
  }
  */

  return;
}

/* ----------------------------------------------------------------------
  calculating total exchange propensity
------------------------------------------------------------------------- */

void Appcoros::ballistic_probability(int n)
{
  int i,k,periodicity[3];
  double dij[4],lprd[3];

  periodicity[0] = domain->xperiodic;
  periodicity[1] = domain->yperiodic;
  periodicity[2] = domain->zperiodic;
  lprd[0] = domain->xprd;
  lprd[1] = domain->yprd;
  lprd[2] = domain->zprd;

  pn_local[n] = pn_global[n] = 0.0;
  for(i = 0; i < nlocal; i ++) {pmix[n][i] = 0.0;}
  for(i = 0; i < nlocal; i ++) {

     for(k = 0; k < 3; k++) {
        dij[k] = xyz[i][k] - xmix[n][k];
        if(periodicity[k] && dij[k] >= lprd[k]/2.0) dij[k] -= lprd[k];
        if(periodicity[k] && dij[k] <= -lprd[k]/2.0) dij[k] += lprd[k];
     }

     dij[3] = sqrt(dij[0]*dij[0]+dij[1]*dij[1]+dij[2]*dij[2]);
     double pb = 0.0;
     if(dij[3] > 0.0) pb = exp(-fabs(dij[3])/rdamp[n]);

     pn_local[n] += pb;
     pmix[n][i] = pn_local[n];
  }

  MPI_Allreduce(&pn_local[n],&pn_global[n],1,MPI_DOUBLE,MPI_SUM,world);
  if(pn_global[n] == 0.0) {
    error->all(FLERR,"total exchange probability returns 0!");
  }
}

/* ----------------------------------------------------------------------
  calculating total energy
------------------------------------------------------------------------- */

double Appcoros::total_energy( )
{
  double penergy = 0.0;
  for(int i = 0; i < nlocal; i++) {

  penergy += sites_energy(i,engstyle);

  }
/*
  if(elastic_flag) {
    for(int j = 0; j < nlocal; j++) {
      int jtype = element[j];
      penergy += elastic_energy(j,jtype);
    }
  }
*/




  return penergy/2.0;
}

/* ----------------------------------------------------------------------
  calculating "average" concentration at each lattice site
  for 1NN interaction (engstyle==1): Ci(a)=0.5*ni(a)+0.5*sum(nj(a))/sum(nj)
  for 2NN interaction (engstyle==2): Ci(a)=0.25*ni(a)+0.5*sum(nj(a))/sum(nj)+0.25*sum(nk(a))/sum(nk))
  subscript i for corners atoms and j for 1NN and k for 2NN of i
  The contributions of corner and center atoms (of bcc cell) are equal
------------------------------------------------------------------------- */

// void Appcoros::concentration_field( )
// {
//   int i,j,jd;
//
//   double *ci = new double[nelement];
//   for(i = 0; i < nlocal; i++) {
//     for(j = 0; j < nelement; j++) ci[j] = 0.0;
//
//     ci[element[i]-1] += 0.25;
//     if(engstyle == 1) ci[element[i]-1] += 0.25;
//
//     for(j = 0; j < numneigh[i]; j ++) {
//       jd = neighbor[i][j];
//       ci[element[jd]-1] += 0.5*1.0/numneigh[i];
//     }
//
//     if(engstyle == 2) {
//       for(j = 0; j < numneigh2[i]; j ++) {
//         jd = neighbor2[i][j];
//         ci[element[jd]-1] += 0.25*1.0/numneigh2[i];
//       }
//     }
//
//     for(j = ndiffusion; j < ndiffusion + concentrationflag; j++) { // count for concentrationflag elements
//       disp[j][i] = ci[j-ndiffusion];
//     }
//   }
//
//   delete [] ci;
// }

/* ----------------------------------------------------------------------
  check if the vacancy is trapped by solute by counting its
  first and second NNs, return 0 if any of those is non-Fe
------------------------------------------------------------------------- */

int Appcoros::vacancy_trap(int i)
{
  int j,jd;

  //energy from 1NN bonds
  for (j = 0; j < numneigh[i]; j++) {
    jd = neighbor[i][j];
    if(element[jd] > 2) return 0;
  }

  //energy from 2NN bonds
  if (engstyle == 2) {
    for (j = 0; j < numneigh2[i]; j++) {
    jd = neighbor2[i][j];
    if(element[jd] > 2) return 0;
    }
  }

  return 1;
}

/* ----------------------------------------------------------------------
  record KMC time and real time
------------------------------------------------------------------------- */
void Appcoros::time_tracer(double dt)
{
  dt_akmc += dt;
  dt_real += dt*itrap; // not contribute to real time if trapped by solute
}

/* ----------------------------------------------------------------------
  check if the vacancy is trapped by solute, contribute 0 to physical
  time if so, return zero if itrap == 1; 1 otherwise
------------------------------------------------------------------------- */

double Appcoros::real_time(double t)
{
  double treal_all,takmc_all;
  double dtreal_all,dtakmc_all;
  double nprcs = (double)(domain->nprocs);

  dtreal_all = dtakmc_all = 0.0;
  MPI_Allreduce(&dt_real,&dtreal_all,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&dt_akmc,&dtakmc_all,1,MPI_DOUBLE,MPI_SUM,world);

  //compute fvt every dt_interval
  itime_current = static_cast<int> (t/dt_interval);
  treal_me += dt_real;
  takmc_me += dt_akmc;

  if(itime_current > itime_old)  {
    treal_all = takmc_all = 0.0;
    MPI_Allreduce(&treal_me,&treal_all,1,MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(&takmc_me,&takmc_all,1,MPI_DOUBLE,MPI_SUM,world);
    if(takmc_all > 0.0) fvt = treal_all/takmc_all/nprcs;
    itime_old = itime_current;
    treal_me = 0.0;
    takmc_me = 0.0;
    dt_real = 0.0;
    dt_akmc = 0.0;
  }

  if(dtakmc_all == 0.0) return 1.0;
  return dtreal_all/dtakmc_all; //averaged real time over all processors
}

/* ----------------------------------------------------------------------
  map n by n matrix to a vector
------------------------------------------------------------------------- */

int Appcoros::ibonde(int a, int b, int c)
{
  return ((a-1)*c + b - a*(a-1)/2);
}

/* ----------------------------------------------------------------------
  grow memory for dislocation
------------------------------------------------------------------------- */
/*
void Appcoros::grow_dislocations()
{
  int n = ndislocation + 1;
  memory->grow(dislocation_type,n,"app/coros:dislocation_type");
  memory->grow(burgers,n,3,"app/coros:burgers");
  memory->grow(xdislocation,n,3,"app/coros:xdislocation");
  memory->grow(line_vector,n,"app/coros:line_vector");
  memory->grow(dislocation_radius,n,"app/coros:dislocation_radius");
  memory->grow(nsegment,n,"app/coros:nsegment");
}
*/
/* ----------------------------------------------------------------------
  grow memory for sink
------------------------------------------------------------------------- */

void Appcoros::grow_sinks()
{
  int n = nsink + 1;
  // memory->grow(sink_shape,n,"app/coros:sink_shape");
  // memory->grow(sink_strength,n,"app/coros:sink_strength");
  // memory->grow(xsink,n,3,"app/coros:xsink");
  // memory->grow(sink_type,n,"app/coros:sink_type");
  // memory->grow(sink_normal,n,"app/coros:sink_normal");
  // memory->grow(sink_segment,n,"app/coros:sink_segmant");
  // memory->grow(sink_radius,n,"app/coros:sink_radius");
  // memory->grow(sink_mfp,n,"app/coros:sink_mfp");
  // memory->grow(nabsorption,n,"app/coros:nabsorption");
  memory->grow(sink_shape,n,"app/ris:sink_shape");
  memory->grow(sink_range,n,"app/ris:sink_range");
  memory->grow(xsink,n,3,"app/ris:xsink");
  memory->grow(sink_normal,n,"app/ris:sink_normal");
  memory->grow(sink_segment,n,"app/ris:sink_segmant");
  memory->grow(sink_radius,n,"app/ris:sink_radius");

}

/* ----------------------------------------------------------------------
  grow memory for reaction
------------------------------------------------------------------------- */
///*
void Appcoros::grow_reactions()
{
  int n = nreaction + 1;
  memory->grow(rsite,n,"app/coros:rsite");
  memory->grow(rinput,n,"app/coros:rinput");
  memory->grow(routput,n,"app/coros:routput");
  memory->grow(rcount,n,"app/coros:rcount");
  memory->grow(renable,n,"app/coros:renable");
  memory->grow(rtarget,n,"app/coros:rtarget");
  memory->grow(rbarrier,n,"app/coros:rbarrier");
  memory->grow(rrate,n,"app/coros:rrate");
  memory->grow(target_local,n,"app/coros:target_local");
  memory->grow(target_global,n,"app/coros:target_global");
}
//*/
/* ----------------------------------------------------------------------
  grow memory for ballistic mixing
------------------------------------------------------------------------- */
/*
void Appcoros::grow_ballistic()
{
  int n = nballistic + 1;
  int m = 3;
  memory->grow(bfreq,n,"app/coros:bfreq");
  memory->grow(time_old,n,"app/coros:time_old");
  memory->grow(time_new,n,"app/coros:time_new");
  memory->grow(rdamp,n,"app/coros:rdamp");
  memory->grow(pn_local,n,"app/coros:pn_local");
  memory->grow(pn_global,n,"app/coros:pn_global");
  memory->grow(xmix,n,m,"app/coros:xmix");
  memory->grow(pmix,n,nlocal,"app/coros:pmix");
}
*/
/* ----------------------------------------------------------------------
  create sinks
------------------------------------------------------------------------- */

void Appcoros::sink_creation(int n)
{
  int i,j,periodicity[3];
  int shape = sink_shape[n];
  int normal = sink_normal[n];
  //int ntype = sink_type[n];
  int segment = sink_segment[n];
  int nlattice = nlocal + nghost;
  double dx,dij[3],rik,rjk,lprd[3];
  double radius = sink_radius[n];
  //double strength = sink_strength[n]*sink_strength[n];

  double range = sink_range[n]*sink_range[n];


  // get periodicity and box length
  periodicity[0] = domain->xperiodic;
  periodicity[1] = domain->yperiodic;
  periodicity[2] = domain->zperiodic;
  lprd[0] = domain->xprd;
  lprd[1] = domain->yprd;
  lprd[2] = domain->zprd;
  // shape =1: straight line sink, e.g., dislocations
  if(shape == 1) {
    for(i = 0; i < nlattice; i++){
      for(j = 0; j < 3; j ++) {
        if(j != normal) dij[j] = xyz[i][j]-xsink[n][j];
        else dij[j] = 0.0;
        if(periodicity[j] && dij[j] >= lprd[j]/2.0) dij[j] -= lprd[j];
        if(periodicity[j] && dij[j] <= -lprd[j]/2.0) dij[j] += lprd[j];
      }
      dx = dij[0]*dij[0]+dij[1]*dij[1]+dij[2]*dij[2];
      //if(dx < strength) isink[i][ntype-1] = 1;
      if(dx < range) isink[i] = n;
    }
  }
  // circular or polygon sinks, e.g., dislocation loops
  else if (shape == 2) {
    if( segment == 0) {
      for(i = 0; i < nlattice; i++){
        dx = 0.0;
        for( j = 0; j < 3; j ++) {
          dij[j] = xyz[i][j]-xsink[n][j];
          if(periodicity[j] && dij[j] >= lprd[j]/2.0) dij[j] -= lprd[j];
          if(periodicity[j] && dij[j] <= -lprd[j]/2.0) dij[j] += lprd[j];
          if(j != normal) dx += dij[j]*dij[j];
        }
        rik = sqrt(dx);
        rjk = (rik-radius)*(rik-radius) + dij[normal]*dij[normal];
        //if( rjk < strength) isink[i][ntype-1] = 1; // LC comment
        if( rjk < range) isink[i] = n;
      }
    }
    else {
      if( segment < 4) error->all(FLERR,"wrong number of segment for polygon sinks");
      double pi = acos(-1.0);
      double theta = 2*pi/segment;
      double theta_xyz;
      int b = normal + 1;
      int c = b + 1;
      if( b > 2) b -= 3;
      if( c > 2) c -= 3;
      for(i = 0; i < nlattice; i++){
        dx = 0.0;
        for( j = 0; j < 3; j ++) {
          dij[j] = xyz[i][j]-xsink[n][j];
          if(periodicity[j] && dij[j] >= lprd[j]/2.0) dij[j] -= lprd[j];
          if(periodicity[j] && dij[j] <= -lprd[j]/2.0) dij[j] += lprd[j];
          if(j != normal) dx += dij[j]*dij[j];
        }
        if(xyz[i][c] == 0 && xyz[i][b] >= 0) theta_xyz = 0.0;
        else if(xyz[i][c] == 0 && xyz[i][b] < 0) theta_xyz = pi;
        else if(xyz[i][b] == 0) theta_xyz = pi/2;
        else theta_xyz = atan(xyz[i][c] / xyz[i][b]);
        if ( theta_xyz < 0 ) theta_xyz += pi;
        int d = static_cast<int>(theta_xyz/theta);
        double phi = fabs(theta_xyz - (d+0.5)*theta);
        rik = sqrt(dx);
        rjk = (rik-radius/cos(phi))*(rik-radius/cos(phi)) + dij[normal]*dij[normal];
        //if( rjk < strength) isink[i][ntype-1] = 1; // LC
        if( rjk < range) isink[i] = n;
      }
    }
  }
  // planar sinks, e.g., grain boundaries
  else if (shape == 3) {
    for(i = 0; i < nlattice; i++){
      dx = xyz[i][normal]-xsink[n][normal];
      if(periodicity[normal] && dx >= lprd[normal]/2.0) dx -= lprd[normal];
      if(periodicity[normal] && dx <= -lprd[normal]/2.0) dx += lprd[normal];
      //if(dx*dx < strength) isink[i][ntype-1] = 1; // LC
      if(dx*dx < range) isink[i] = n;
    }
  }
  // 3D spherical sinks, e.g., precipitates
  else {
    for(i = 0; i < nlattice; i++){
      for( j = 0; j < 3; j ++) {
        dij[j] = xyz[i][j]-xsink[n][j];
        if(periodicity[j] && dij[j] >= lprd[j]/2.0) dij[j] -= lprd[j];
        if(periodicity[j] && dij[j] <= -lprd[j]/2.0) dij[j] += lprd[j];
      }
      dx = dij[0]*dij[0]+dij[1]*dij[1]+dij[2]*dij[2];
      //if(dx < strength) isink[i][ntype-1] = 1; // LC
      if(dx < range) isink[i] = n;
    }
  }
}
/* ----------------------------------------------------------------------
  sink_statistics, record the # sites and their ids for each sink
------------------------------------------------------------------------- */
//
// void Appcoros::sink_statistics()
// { int i,j,maxsite,id;
//
//   // Record the sink site ids
//
//   maxsite = 0;
//   memory->create(sinksite,nsink,"app/seg:sinksite"); // # of sink sites for each sink
//   for(i = 0; i < nsink; i++) {
//      sinksite[i] = 0;
//      for(j = 0; j < nlocal; j++) {
//         if(isink[j] == i+1) sinksite[i] ++; // sink id starts from 1
//      }
//      if(sinksite[i] > maxsite) maxsite = sinksite[i];
//   }
//
//   memory->create(sinkid,nsink,maxsite,"app/seg:sinkid"); // ids of each sink site for each sink
//   for(i = 0; i < nsink; i++) {
//      id = 0;
//      for(j = 0; j < nlocal; j++) {
//         if(isink[j] == i+1) {
//           sinkid[i][id] = j;
// 	  id++;
// 	}
//      }
//   }
//
//   return;
// }
/* ----------------------------------------------------------------------
  calculate dislocation stess field
------------------------------------------------------------------------- */
/*
void Appcoros::stress_field(int n)
{
  if(dislocation_type[n] == 1) stress_dislocation(n); //straight dislocation
  else stress_loop(n); //loop
}
*/
/* ----------------------------------------------------------------------
  calculate stess field for straight dislocations
------------------------------------------------------------------------- */
/*
void Appcoros::stress_dislocation(int n)
{ int i,j,k,l,ii;
  int nlattice = nlocal + nghost;
  int periodicity[3];
  double pix;
  double ni[3],m[3],nn[3][3],nm[3][3],nni[3][3];
  double temp1[3][3],temp2[3][3];
  double bmatx[3][3],qmatx[3][3],smatx[3][3];
  double sigma[3][3],stran[3][3],stres[3][3];
  double tvect[3],dij[3];
  double lprd[3];
  // get periodicity and box length
  periodicity[0] = domain->xperiodic;
  periodicity[1] = domain->yperiodic;
  periodicity[2] = domain->zperiodic;
  lprd[0] = domain->xprd;
  lprd[1] = domain->yprd;
  lprd[2] = domain->zprd;
  // calculate the elastic tensor
  elastic_tensor();
  for (i = 0; i < 3; i ++) {
    if(i == line_vector[n]) tvect[i] = 1.0;
    else tvect[i] = 0.0;
  }
  pix = acos(-1.0);
  vector_normalize(tvect);
  stroh(tvect,qmatx,bmatx,smatx);
  // calculate stress at each lattice site
  for( int natom = 0; natom < nlattice; natom ++) {
    for( j = 0; j < 3; j ++) {
      dij[j] = xyz[natom][j] - xdislocation[n][j];
      if(periodicity[j] && dij[j] > lprd[j]/2.0) dij[j] -= lprd[j];
      if(periodicity[j] && dij[j] <= -lprd[j]/2.0) dij[j] += lprd[j];
    }
    double norm = dij[0]*tvect[0] + dij[1]*tvect[1] + dij[2]*tvect[2];
    for(i = 0; i < 3; i ++) m[i] = dij[i] - tvect[i]*norm;
    double d = sqrt(m[0]*m[0] + m[1]*m[1] + m[2]*m[2]);
    if(d == 0.0) {
      for(i = 0; i < 6; i ++) stress[natom][i] = 0.0;
    } else {
      vector_normalize(m);
      cross_product(tvect,m,ni);
      for(j = 0; j < 3; j ++) {
        for(k = 0; k < 3; k ++) {
          nn[j][k] = 0.0;
          nm[j][k] = 0.0;
        }
      }
      for(i = 0; i < 3; i ++) {
        for(j = 0; j < 3; j ++) {
          for(k = 0; k < 3; k ++) {
            for(l = 0; l < 3; l ++) {
              nn[i][j] += ni[k]*cijkl[i][j][k][l]*ni[l];
              nm[i][j] += ni[k]*cijkl[i][j][k][l]*m[l];
            }
          }
        }
      }
      matrix_inversion(nn,nni);
      // inteimediate matrices
      for(i = 0; i < 3; i ++) {
        for(j = 0; j < 3; j ++) {
          temp1[i][j] = bmatx[i][j];
          for(k = 0; k < 3; k ++) {
            temp1[i][j] += nm[i][k]*smatx[k][j];
          }
       }
      }
      for(i = 0; i < 3; i ++) {
        for(j = 0; j < 3; j ++) {
          temp2[i][j] = 0.0;
          for(k = 0; k < 3; k ++) {
            temp2[i][j] += nni[i][k]*temp1[k][j];
          }
        }
      }
      // form sigma angular
      for(i = 0; i < 3; i ++) {
        for(j = 0; j < 3; j ++) {
          sigma[i][j] = 0.0;
          stran[i][j] = 0.0;
          stres[i][j] = 0.0;
        }
      }
      for(i = 0; i < 3; i ++) {
        for(j = 0; j < 3; j ++) {
          for(k = 0; k < 3; k ++) {
            for(l = 0; l < 3; l ++) {
              for(ii = 0; ii < 3; ii ++) {
                sigma[i][j] += cijkl[i][j][k][l]*burgers[n][ii]*
                   (-m[l]*smatx[k][ii] + ni[l]*temp2[k][ii])/2/pix;
              }
            }
          }
        }
      }
      for(i = 0; i < 3; i ++) {
        for(j = 0; j < 3; j ++) {
          // stres[i][j] = sigma[i][j]/d;
          stres[i][j] = sigma[i][j]*d/(d*d + dcore*dcore); //stress converge at core distance
        }
      }
    // calculate strain and stress
      for(i = 0; i < 3; i ++) {
        for(j = 0; j < 3; j ++) {
          for(ii = 0; ii < 3; ii ++) {
            // stran[i][j] += burgers[n][ii]*(-m[j]*smatx[i][ii] + ni[j]*temp2[i][ii])/d;
            stran[i][j] += burgers[n][ii]*(-m[j]*smatx[i][ii] + ni[j]*temp2[i][ii])*d/(d*d + dcore*dcore);
          }
        }
      }
      for(i = 0; i < 3; i ++) {
        for(j = 0; j < 3; j ++) {
          stres[i][j] = 0.0;
          for(k = 0; k < 3; k ++) {
            for(l = 0; l < 3; l ++) {
              stres[i][j] += cijkl[i][j][k][l]*stran[k][l];
            }
          }
        }
      }
      // record stress at lattice point
      stress[natom][0] = stres[0][0];
      stress[natom][1] = stres[1][1];
      stress[natom][2] = stres[2][2];
      stress[natom][3] = (stres[0][1]+stres[1][0])/2.0;
      stress[natom][4] = (stres[0][2]+stres[2][0])/2.0;
      stress[natom][5] = (stres[1][2]+stres[2][1])/2.0;
    }
  }
}
*/
/* ----------------------------------------------------------------------
  calculate stess field for dislocation loops
------------------------------------------------------------------------- */
/*
void Appcoros::stress_loop(int n)
{
  int i,j;
  int normal,nseg,iseg,jseg;
  int nlattice = nlocal + nghost;
  int periodicity[3];
  double pix,theta;
  double dij[3],xseg[40][3]; //maximum 40 segment to represent a circle
  double A[3],B[3],P[3];
  double bv[3],rloop;  // burgers vector and loop normal & radius
  double stres[3][3],sstres[3][3]; // stress tensor at lattice site
  double lprd[3];
  // get periodicity and box length
  periodicity[0] = domain->xperiodic;
  periodicity[1] = domain->yperiodic;
  periodicity[2] = domain->zperiodic;
  lprd[0] = domain->xprd;
  lprd[1] = domain->yprd;
  lprd[2] = domain->zprd;
  elastic_tensor();
  for(i = 0; i < 3; i ++) {
    bv[i] = burgers[n][i];
  }
  rloop = dislocation_radius[n];
  normal = line_vector[n];
  nseg = nsegment[n];
  pix = acos(-1.0);
  theta = 2*pix/nseg;
  int b = normal + 1;
  int c = b + 1;
  if( b > 2) b -= 3;
  if( c > 2) c -= 3;
  // get vertice at intersections of segments
  for (iseg = 0; iseg < nseg+1; iseg ++) {
    xseg[iseg][b] = rloop*cos(theta*iseg);
    xseg[iseg][c] = rloop*sin(theta*iseg);
    xseg[iseg][normal] = 0.0;
  }
  // calculate stress at each lattice site
  for( int natom = 0; natom < nlattice; natom ++) {
    for( j = 0; j < 3; j ++) {
      dij[j] = xyz[natom][j] - xdislocation[n][j];
      if(periodicity[j] && dij[j] >= lprd[j]/2.0) dij[j] -= lprd[j];
      if(periodicity[j] && dij[j] <= -lprd[j]/2.0) dij[j] += lprd[j];
    }
    for( i = 0; i < 3; i++) {
      P[i] = dij[i];
      for( j = 0; j < 3; j ++) stres[i][j] = 0.0;
    }
    for( iseg = 0; iseg < nseg; iseg ++) {
      jseg = iseg+1;
      if(jseg == nseg) jseg = 0;
      for(i = 0; i < 3; i ++) {
        A[i] = xseg[iseg][i];
        B[i] = xseg[jseg][i];
        for( j = 0; j < 3; j ++) sstres[i][j] = 0.0;
      }
      seg_stress(A,B,P,bv,sstres);
      for(i = 0; i < 3; i ++) {
        for(j = 0; j < 3; j ++) stres[i][j] += sstres[i][j];
      }
    }
    // record stress at lattice point
    stress[natom][0] = stres[0][0];
    stress[natom][1] = stres[1][1];
    stress[natom][2] = stres[2][2];
    stress[natom][3] = (stres[0][1]+stres[1][0])/2.0;
    stress[natom][4] = (stres[0][2]+stres[2][0])/2.0;
    stress[natom][5] = (stres[1][2]+stres[2][1])/2.0;
  }
}
*/
/* ----------------------------------------------------------------------
  Calculate stress field due to segment AB
------------------------------------------------------------------------- */
/*
void Appcoros::seg_stress( double A[3], double B[3], double P[3], double bv[3], double sstres[3][3])
{ int i,j;
  double xi[3],PA[3],PB[3],x[3],ni[3];
  double norm,eps,dotx,beta1,beta2;
  double tau1[3],tau2[3],m1[3],m2[3];
  double sigma1[3][3],sigma2[3][3],sigma_p1[3][3],sigma_p2[3][3];
  eps = 1.0e-6;
  for(i = 0; i < 3; i ++) {
    x[i] = 0.0;
    xi[i] = B[i] - A[i];
  }
  vector_normalize(xi);
  for(i = 0; i < 3; i ++) {
    PA[i] = P[i] - A[i];
    PB[i] = P[i] - B[i];
  }
  dotx = 0.0;
  for(i = 0; i < 3; i ++) dotx += PA[i]*xi[i];
  for(i = 0; i < 3; i ++) x[i] += PA[i] - dotx*xi[i];
  double normx = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
  // return zero if co-linear
  if(normx < eps) {
    for(i = 0; i < 3; i ++) {
      for(j = 0; j < 3; j ++) sstres[i][j] = 0.0;
    }
    return;
    // if normx < dcore, scale normx to dcore
  } else if (normx < dcore) {
    for(i = 0; i < 3; i ++) {
      PA[i] = dotx*xi[i] + x[i]*dcore/normx;
      PB[i] = PA[i] + A[i] - B[i];
    }
  }
  // ABP plane normal to ni
  vector_normalize(x);
  norm = sqrt(PA[0]*PA[0] + PA[1]*PA[1] + PA[2]*PA[2]);
  for(i = 0; i < 3; i ++) tau1[i] = PA[i]/norm;
  norm = sqrt(PB[0]*PB[0] + PB[1]*PB[1] + PB[2]*PB[2]);
  for(i = 0; i < 3; i ++) tau2[i] = PB[i]/norm;
  // ni & m1 & m2 vactor
  cross_product(xi,x,ni);
  cross_product(ni,tau1,m1);
  cross_product(ni,tau2,m2);
  // angles
  dotx = 0.0;
  for(i = 0; i < 3; i ++)  dotx += tau1[i]*xi[i];
  beta1 = acos(dotx);
  dotx = 0.0;
  for(i = 0; i < 3; i ++)  dotx += tau2[i]*xi[i];
  beta2 = acos(dotx);
  // Augular factors & stress factors
  sigma_A(tau1, m1, bv, sigma1);
  sigma_A(tau2, m2, bv, sigma2);
  sigma_P(tau1, ni, bv, sigma_p1);
  sigma_P(tau2, ni, bv, sigma_p2);
  // calculate the stress field due to AB segment
  for(i = 0; i < 3; i ++) {
    for(j = 0; j < 3; j ++) {
      sstres[i][j] = (cos(beta1)*sigma1[i][j] - cos(beta2)*sigma2[i][j])*normx/2
                     /(normx*normx + dcore*dcore) + (sin(beta2)*sigma_p2[i][j] - sin(beta1)*
                     sigma_p1[i][j])*normx/2/(normx*normx + dcore*dcore);
    }
  }
}
*/
/* ----------------------------------------------------------------------
  Calculate sigma_A, used in seg_stress
------------------------------------------------------------------------- */
/*
void Appcoros::sigma_A(double t[3], double m[3], double bv[3], double sigma[3][3])
{
  int i,j,k,l,ii;
  double pix;
  double ni[3],nx[3],tx[3];
  double bmatx[3][3],qmatx[3][3],smatx[3][3];
  double bpmat[3][3],qpmat[3][3],spmat[3][3];
  double nn[3][3],nm[3][3],nni[3][3];
  double temp[3][3],ma[3][3];
  pix = acos(-1.0);
  vector_normalize(t);
  cross_product(t,m,ni);
  vector_normalize(ni);
  // qpmatx,bpmatx,spmatx
  for(i = 0; i < 3; i ++) {
    nx[i] = ni[i];
    tx[i] = t[i];
  }
  stroh_p(tx,nx,qpmat,bpmat,spmat);
  // qmatx,bmatx,smatx
  for(i = 0; i < 3; i ++) tx[i] = t[i];
  stroh(tx,qmatx,bmatx,smatx);
  // nn & nm
  for(i = 0; i < 3; i ++) {
    for(j = 0; j < 3; j ++) {
      nn[i][j] = 0.0;
      nm[i][j] = 0.0;
    }
  }
  for(j = 0; j < 3; j ++) {
    for(k = 0; k < 3; k ++) {
      for(i = 0; i < 3; i ++) {
        for(l = 0; l < 3; l ++) {
          nn[j][k] += ni[i]*cijkl[i][j][k][l]*ni[l];
          nm[j][k] += ni[i]*cijkl[i][j][k][l]*m[l];
        }
      }
    }
  }
  matrix_inversion(nn,nni);
  //temporaty matrix ma
  for(i = 0; i < 3; i ++) {
    for(j = 0; j < 3; j ++) {
      temp[i][j] = 0.0 ;
      for(k = 0; k < 3; k ++)  temp[i][j] += nm[i][k]*smatx[k][j];
    }
  }
  for(i = 0; i < 3; i ++) {
    for(j = 0; j < 3; j ++)  temp[i][j] += bmatx[i][j];
  }
  for(i = 0; i < 3; i ++) {
    for(j = 0; j < 3; j ++) {
      ma[i][j] = 0.0 ;
      for(k = 0; k < 3; k ++)  ma[i][j] += nni[i][k]*temp[k][j];
    }
  }
  // calculate the angular stress factors
  for(i = 0; i < 3; i ++) {
    for(j = 0; j < 3; j ++)  sigma[i][j] = 0.0;
  }
  for(i = 0; i < 3; i ++) {
    for(j = 0; j < 3; j ++) {
      for(k = 0; k < 3; k ++) {
        for(l = 0; l < 3; l ++) {
          for(ii = 0; ii < 3; ii ++) {
            sigma[i][j] += cijkl[i][j][k][l]*bv[ii]*
                 (-m[l]*smatx[k][ii] + ni[l]*ma[k][ii])/2/pix;
          }
        }
      }
    }
  }
}
*/
/* ----------------------------------------------------------------------
  Calculate sigma_P, used in seg_stress
------------------------------------------------------------------------- */
/*
void Appcoros::sigma_P(double t[3], double ni[3], double bv[3], double sigmap[3][3])
{ int i, j, k, l, ii;
  double pix;
  double nn[3][3],nm[3][3],nt[3][3],nni[3][3],mi[3];
  double bmatx[3][3],qmatx[3][3],smatx[3][3];
  double bpmat[3][3],qpmat[3][3],spmat[3][3];
  double temp1[3][3],temp2[3][3],temp3[3][3];
  double ma[3][3];
  pix = acos(-1.0);
  vector_normalize(t);
  vector_normalize(ni);
  cross_product(ni,t,mi);
  stroh_p(t,ni,qpmat,bpmat,spmat);
  stroh(t,qmatx,bmatx,smatx);
  // nn & nm
  for(i = 0; i < 3; i ++) {
    for(j = 0; j < 3; j ++) {
      nn[i][j] = 0.0;
      nm[i][j] = 0.0;
      nt[i][j] = 0.0;
    }
  }
  for(i = 0; i < 3; i ++) {
    for(j = 0; j < 3; j ++) {
      for(k = 0; k < 3; k ++) {
        for(l = 0; l < 3; l ++) {
          nn[j][k] += ni[i]*cijkl[i][j][k][l]*ni[l];
          nm[j][k] += ni[i]*cijkl[i][j][k][l]*mi[l];
          nt[j][k] += ni[i]*cijkl[i][j][k][l]*t[l];
        }
      }
    }
  }
  matrix_inversion(nn,nni);
  // temporaty matrix ma
  for(i = 0; i < 3; i ++) {
    for(j = 0; j < 3; j ++) {
      temp1[i][j] = 0.0 ;
      for(k = 0; k < 3; k ++)  temp1[i][j] += nt[i][k]*smatx[k][j];
    }
  }
  for(i = 0; i < 3; i ++) {
    for(j = 0; j < 3; j ++) {
      temp2[i][j] = 0.0 ;
      for(k = 0; k < 3; k ++)  temp2[i][j] += nm[i][k]*spmat[k][j];
    }
  }
  for(i = 0; i < 3; i ++) {
    for(j = 0; j < 3; j ++)    temp3[i][j] = bpmat[i][j] + temp2[i][j] - temp1[i][j];
  }
  for(i = 0; i < 3; i ++) {
    for(j = 0; j < 3; j ++) {
      ma[i][j] = 0.0 ;
      for(k = 0; k < 3; k ++)  ma[i][j] += nni[i][k]*temp3[k][j];
    }
  }
  // calculate the angular stress factors
  for(i = 0; i < 3; i ++) {
    for(j = 0; j < 3; j ++)    sigmap[i][j] = 0.0;
  }
  for(i = 0; i < 3; i ++) {
    for(j = 0; j < 3; j ++) {
      for(k = 0; k < 3; k ++) {
        for(l = 0; l < 3; l ++) {
          for(ii = 0; ii < 3; ii ++) {
            sigmap[i][j] += cijkl[i][j][k][ii]*bv[l]*
                 (t[ii]*smatx[k][l] - mi[ii]*spmat[k][l] + ni[ii]*ma[k][l])/2/pix;
          }
        }
      }
    }
  }
}
*/
/* ----------------------------------------------------------------------
  stress calculation based on stroh formulism, used in stress_dislocation and sigma_A and sigma_P
------------------------------------------------------------------------- */
/*
void Appcoros::stroh( double tvect[3], double qmatx[3][3], double bmatx[3][3], double smatx[3][3])
{
  int i,j,k,l;
  double pix,omega,domega;
  double mvect[3],nvect[3],mmvect[3],nnvect[3];
  double nn[3][3],mm[3][3],nm[3][3],mn[3][3];
  double nni[3][3],nn2[3][3],nn3[3][3];
  pix = acos(-1.0);
  domega = 2*pix/ninteg;
  right_hand_coord(mvect,nvect,tvect);
  for(i = 0; i < 3; i ++) {
    for(j = 0; j < 3; j ++) {
      qmatx[i][j] = 0.0;
      smatx[i][j] = 0.0;
      bmatx[i][j] = 0.0;
    }
  }
  for(int integ = 0; integ < ninteg; integ ++) {
    for(i = 0; i < 3; i ++) {
      for(j = 0; j < 3; j ++) {
        nn[i][j] = 0.0;
        nm[i][j] = 0.0;
        mn[i][j] = 0.0;
        mm[i][j] = 0.0;
      }
    }
    omega = integ*domega;
    for(i = 0; i < 3; i ++) {
      mmvect[i] = mvect[i]*cos(omega)+nvect[i]*sin(omega);
      nnvect[i] = -mvect[i]*sin(omega)+nvect[i]*cos(omega);
    }
    // different with Bulent's version, double check
    for(i = 0; i < 3; i++) {
      for(j = 0; j < 3; j++) {
        for(k = 0; k < 3; k++) {
          for(l = 0; l < 3; l++) {
             nn[i][j] += nnvect[k]*cijkl[i][j][k][l]*nnvect[l];
             nm[i][j] += nnvect[k]*cijkl[i][j][k][l]*mmvect[l];
             mm[i][j] += mmvect[k]*cijkl[i][j][k][l]*mmvect[l];
          }
        }
      }
    }
    for(i = 0; i < 3; i++) {
      for(j = 0; j < 3; j++) mn[j][i] = nm[i][j];
    }
    matrix_inversion(nn,nni);
    // get nn2
    for(i = 0; i < 3; i++) {
      for(j = 0; j < 3; j++) {
        nn2[i][j] = 0.0;
        for(k = 0; k < 3; k++)  nn2[i][j] += nni[i][k]*nm[k][j];
      }
    }
    // get nn3
    for(i = 0; i < 3; i++) {
      for(j = 0; j < 3; j++) {
        nn3[i][j] = mm[i][j];
        for(k = 0; k < 3; k++)  nn3[i][j] -= mn[i][k]*nn2[k][j];
      }
    }
    // calculate qmatx, smatx and bmatx
    for(i = 0; i < 3; i++) {
      for(j = 0; j < 3; j++) {
        qmatx[i][j] -= nni[i][j]/ninteg;
        smatx[i][j] -= nn2[i][j]/ninteg;
        bmatx[i][j] += nn3[i][j]/ninteg;
      }
    }
  } // loop for integ = 0 - ninteg
}
*/
/* ----------------------------------------------------------------------
  angular factors in stress calculation based on stroh formulism,, used in stress_dislocation and sigma_A and sigma_P
------------------------------------------------------------------------- */
/*
void Appcoros::stroh_p( double t[3], double n0[3], double qpmat[3][3], double bpmat[3][3], double spmat[3][3])
{ int i,j,k,l,ii;
  double pix,value;
  double N1[3],M1[3],ni[3],mi[3];
  double domega,omega,omegas[100];
  double qmatx[3][3],bmatx[3][3],smatx[3][3];
  double nn[3][3],mm[3][3],mn[3][3],nm[3][3],nni[3][3];
  double nt[3][3],tn[3][3],mt[3][3],tm[3][3];
  double F[3][3],Bi[3][3],Si[3][3],Qi[3][3];
  double temp1[3][3],temp2[3][3],temp3[3][3],temp4[3][3];
  double BP[9][100],SP[9][100],QP[9][100];
  pix = acos(-1.0);
  vector_normalize(t);
  vector_normalize(n0);
  for(i = 0; i < 3; i ++) N1[i] = n0[i];
  cross_product(N1,t,M1);
  stroh(t,qmatx,bmatx,smatx);
  // angles
  domega = 2.0*pix/ninteg;
  for(i = 0; i < ninteg+1; i ++)  omegas[i] = i*domega;
  for(i = 0; i < 9; i ++) {
    for(j = 0; j < ninteg+1; j ++) {
      BP[i][j] = 0.0;
      SP[i][j] = 0.0;
      QP[i][j] = 0.0;
    }
  }
  for(i = 0; i < ninteg+1; i ++) {
    omega = omegas[i];
    for(j = 0; j < 3; j ++) {
      mi[j] = M1[j]*cos(omega) + N1[j]*sin(omega);
      ni[j] = -M1[j]*sin(omega) + N1[j]*cos(omega);
    }
    for(ii = 0; ii < 3; ii ++) {
      for(j = 0; j < 3; j ++) {
        nn[ii][j] = 0.0;
        nm[ii][j] = 0.0;
        mm[ii][j] = 0.0;
        mn[ii][j] = 0.0;
        nt[ii][j] = 0.0;
        tn[ii][j] = 0.0;
        mt[ii][j] = 0.0;
        tm[ii][j] = 0.0;
      }
    }
    for(ii = 0; ii < 3; ii ++){
      for(j = 0; j < 3; j ++){
        for(k = 0; k < 3; k ++){
          for(l = 0; l < 3; l ++){
            nn[j][k] += ni[ii]*cijkl[ii][j][k][l]*ni[l];
            nm[j][k] += ni[ii]*cijkl[ii][j][k][l]*mi[l];
            mm[j][k] += mi[ii]*cijkl[ii][j][k][l]*mi[l];
            mn[j][k] += mi[ii]*cijkl[ii][j][k][l]*ni[l];
            nt[j][k] += ni[ii]*cijkl[ii][j][k][l]*t[l];
            tn[j][k] += t[ii]*cijkl[ii][j][k][l]*ni[l];
            mt[j][k] += mi[ii]*cijkl[ii][j][k][l]*t[l];
            tm[j][k] += t[ii]*cijkl[ii][j][k][l]*mi[l];
          }
        }
      }
    }
    matrix_inversion(nn,nni);
    // define F matrix
    for(j = 0; j < 3; j ++) {
      for(k = 0; k < 3; k ++) {
        temp1[j][k] = nt[j][k] + tn[j][k];
      }
    }
    for(j = 0; j < 3; j ++) {
      for(k = 0; k < 3; k ++) {
        temp2[j][k] = 0.0;
        for(l = 0; l < 3; l ++) temp2[j][k] += temp1[j][l]*nni[l][k];
      }
    }
    for(j = 0; j < 3; j ++) {
      for(k = 0; k < 3; k ++) {
        F[j][k] = 0.0;
        for(l = 0; l < 3; l ++) F[j][k] += nni[j][l]*temp2[l][k];
      }
    }
    // define Si matrix
    for(j = 0; j < 3; j ++) {
      for(k = 0; k < 3; k ++) {
        temp1[j][k] = tm[j][k]*sin(omega) - nt[j][k]*cos(omega);
      }
    }
    for(j = 0; j < 3; j ++) {
      for(k = 0; k < 3; k ++) {
        temp2[j][k] = 0.0;
        for(l = 0; l < 3; l ++) temp2[j][k] += nni[j][l]*temp1[l][k];
      }
    }
    for(j = 0; j < 3; j ++) {
      for(k = 0; k < 3; k ++) {
        temp1[j][k] = 0.0;
        for(l = 0; l < 3; l ++) temp1[j][k] += F[j][l]*nm[l][k];
      }
    }
    for(j = 0; j < 3; j ++) {
      for(k = 0; k < 3; k ++) {
        Si[j][k] = -temp1[j][k]*sin(omega) + temp2[j][k];
      }
    }
    // define Qi
    for(j = 0; j < 3; j ++) {
      for(k = 0; k < 3; k ++) {
        Qi[j][k] = -F[j][k]*sin(omega);
      }
    }
    // define Bi
    for(j = 0; j < 3; j ++) {
      for(k = 0; k < 3; k ++) {
        temp1[j][k] = tm[j][k]*sin(omega) - nt[j][k]*cos(omega);
      }
    }
    for(j = 0; j < 3; j ++) {
      for(k = 0; k < 3; k ++) {
        temp2[j][k] = 0.0;
        for(l = 0; l < 3; l ++) temp2[j][k] += nni[j][l]*temp1[l][k];
      }
    }
    // keep this temp1
    for(j = 0; j < 3; j ++) {
      for(k = 0; k < 3; k ++) {
        temp1[j][k] = 0.0;
        for(l = 0; l < 3; l ++) temp1[j][k] += mn[j][l]*temp2[l][k];
      }
    }
    for(j = 0; j < 3; j ++) {
      for(k = 0; k < 3; k ++) {
        temp2[j][k] = 0.0;
        for(l = 0; l < 3; l ++) temp2[j][k] += nni[j][l]*nm[l][k];
      }
    }
    for(j = 0; j < 3; j ++) {
      for(k = 0; k < 3; k ++) {
        temp3[j][k] = tn[j][k]*cos(omega) - mt[j][k]*sin(omega);
      }
    }
    // keep this temp4
    for(j = 0; j < 3; j ++) {
      for(k = 0; k < 3; k ++) {
        temp4[j][k] = 0.0;
        for(l = 0; l < 3; l ++) temp4[j][k] += temp3[j][l]*temp2[l][k];
      }
    }
    for(j = 0; j < 3; j ++) {
      for(k = 0; k < 3; k ++) {
        temp2[j][k] = 0.0;
        for(l = 0; l < 3; l ++) temp2[j][k] += F[j][l]*nm[l][k]*sin(omega);
      }
    }
    // keep this temp3
    for(j = 0; j < 3; j ++) {
      for(k = 0; k < 3; k ++) {
        temp3[j][k] = 0.0;
        for(l = 0; l < 3; l ++) temp3[j][k] += mn[j][l]*temp2[l][k];
      }
    }
    // keep this temp2
    for(j = 0; j < 3; j ++) {
      for(k = 0; k < 3; k ++) {
        temp2[j][k] = mt[j][k] + tm[j][k];
      }
    }
    for(j = 0; j < 3; j ++) {
      for(k = 0; k < 3; k ++) {
        Bi[j][k] = -temp2[j][k]*cos(omega) + temp3[j][k] + temp4[j][k] - temp1[j][k];
      }
    }
    // define BP, SP, QP
    for(j = 0; j < 3; j ++) {
      for(k = 0; k < 3; k ++) {
        BP[3*j+k][i] = Bi[j][k];
        SP[3*j+k][i] = Si[j][k];
        QP[3*j+k][i] = Qi[j][k];
      }
    }
  } // end loop for i = 0 - ninteg
  value = 0.0;
  for(j = 0; j < 3; j ++) {
    for(k = 0; k < 3; k ++) {
      trapozidal(omegas,BP,j,k,value);
      bpmat[j][k] = value/2/pix;
      trapozidal(omegas,SP,j,k,value);
      spmat[j][k] = value/2/pix;
      trapozidal(omegas,QP,j,k,value);
      qpmat[j][k] = value/2/pix;
    }
  }
}
*/
/* ----------------------------------------------------------------------
  establish elastic tensor matrix from c11 c12 & c44, used in stress_loop and stress_dislocation
------------------------------------------------------------------------- */
/*
void Appcoros::elastic_tensor()
{
  int i,j,k,l;
  double anisotropy,delta[3][3];
  // shift c44
  anisotropy = fabs(2*c44 + c12 - c11);
  if(anisotropy < 1.0e-6) c44 += 1.0e-6;
  // identity matrix
  for(i = 0; i < 3; i++) {
    for(j = 0; j < 3; j++) {
      if(i == j) delta[i][j] = 1;
      else delta[i][j] = 0;
    }
  }
  // cijkl tensor
  anisotropy = 2*c44 + c12 - c11;
  for(i = 0; i < 3; i++) {
    for(j = 0; j < 3; j++) {
      for(k = 0; k < 3; k++) {
        for(l = 0; l < 3; l++) {
          cijkl[i][j][k][l] = c44*(delta[i][k]*delta[j][l] + delta[i][l]*delta[j][k])
                              + c12*delta[i][j]*delta[k][l]
                              -anisotropy*delta[i][j]*delta[k][l]*delta[i][k];
        }
      }
    }
  }
}
*/
/* ----------------------------------------------------------------------
  3x3 matrix inversion
------------------------------------------------------------------------- */

void Appcoros::matrix_inversion(double mm[3][3], double nn[3][3])
{
   double det = mm[0][0]*(mm[2][2]*mm[1][1] - mm[2][1]*mm[1][2])
                - mm[1][0]*(mm[2][2]*mm[0][1] - mm[2][1]*mm[0][2])
                + mm[2][0]*(mm[1][2]*mm[0][1] - mm[1][1]*mm[0][2]);

   if(det <= 0) error->all(FLERR,"matrix not invertible!");

   // get the inverse matrix

   nn[0][0] = (mm[2][2]*mm[1][1] - mm[2][1]*mm[1][2])/det;
   nn[0][1] = -(mm[2][2]*mm[0][1] - mm[2][1]*mm[0][2])/det;
   nn[0][2] = (mm[1][2]*mm[0][1] - mm[1][1]*mm[0][2])/det;
   nn[1][0] = -(mm[2][2]*mm[1][0] - mm[2][0]*mm[1][2])/det;
   nn[1][1] = (mm[2][2]*mm[0][0] - mm[2][0]*mm[0][2])/det;
   nn[1][2] = -(mm[1][2]*mm[0][0] - mm[1][0]*mm[0][2])/det;
   nn[2][0] = (mm[2][1]*mm[1][0] - mm[2][0]*mm[1][1])/det;
   nn[2][1] = -(mm[2][1]*mm[0][0] - mm[2][0]*mm[0][1])/det;
   nn[2][2] = (mm[1][1]*mm[0][0] - mm[1][0]*mm[0][1])/det;

}

/* ----------------------------------------------------------------------
  cross product of double vectors
------------------------------------------------------------------------- */

void Appcoros::cross_product( double m[3], double n[3], double l[3])
{
  l[0] = m[1]*n[2] - m[2]*n[1];
  l[1] = m[2]*n[0] - m[0]*n[2];
  l[2] = m[0]*n[1] - m[1]*n[0];

}

/* ----------------------------------------------------------------------
  cross product of int vectors
------------------------------------------------------------------------- */

void Appcoros::cross_product( int m[3], int n[3], int l[3])
{
  l[0] = m[1]*n[2] - m[2]*n[1];
  l[1] = m[2]*n[0] - m[0]*n[2];
  l[2] = m[0]*n[1] - m[1]*n[0];

}

/* ----------------------------------------------------------------------
  normalize double vector
------------------------------------------------------------------------- */

void Appcoros::vector_normalize( double m[3])
{
  double norm = sqrt(m[0]*m[0] + m[1]*m[1] + m[2]*m[2]);
  if(norm == 0) error->all(FLERR,"can not normalize zero vector");
  for(int i = 0; i < 3; i ++) m[i] /= norm;

}

/* ----------------------------------------------------------------------
  right-hand coordination system, m x n x t
------------------------------------------------------------------------- */

void Appcoros::right_hand_coord( double m[3], double n[3], double t[3])
{
  double taux[3];
  // random vector taux to find m
  taux[0] = rancoros->uniform();
  taux[1] = rancoros->uniform();
  taux[2] = rancoros->uniform();

  vector_normalize(t);
  vector_normalize(taux);

  // m = t x taux
  cross_product(t,taux,m);
  vector_normalize(m);

  // n = t x m
  cross_product(t,m,n);
  vector_normalize(n);

}

/* ----------------------------------------------------------------------
  trapozidal integration
------------------------------------------------------------------------- */

void Appcoros::trapozidal(double omegas[100], double XP[9][100], int j, int k, double value)
{ int i, ii;
  double dy,domega;

  for(i = 0; i < ninteg; i ++) {
    ii = i + 1;
    domega = omegas[ii] - omegas[i];
    dy = XP[3*j+k][ii] + XP[3*j+k][i];
    value += dy*domega/2.0;
  }
}

/* ----------------------------------------------------------------------
  Below functions are created by LC
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
  update region; comment but keep
------------------------------------------------------------------------- */
/*
void Appcoros::update_region(int i,int j, int r)
{
  //define varialbes
  n_update_list = 0;
  //update_list  // this is the list we want to check, already defined in .h file
  //check type/element and update type
  if(r == 1 || r == 3 || r == 4){
    //if(type[i] == 1 && type[j] == 1){}//diffusion in bulk, no need to update
    //diffusion from bulk to surface //need to check neighbor
    if(element[i] != 2 && element[j] == 2 && type[i] == 1 && type[j] == 2){
      type[i] = 2;
      type[j] = 3;
      update_list[0] = j;
      update_list[1] = i;
      n_update_list = 2;
    }
    // for surface diffusion
    if(element[i] == 2 && element[j] != 2 && type[i] == 2 && type[j] == 3){
      type[i] = 3;
      type[j] = 2;
      update_list[0] = j;
      update_list[1] = i;
      n_update_list = 2;
      while (n_update_list != 0){
      n_update_list = update_surface_diff(n_update_list); //!!check if need to combine with update_neighbor_check function
    }
    }
  }
  if(r == 2){//reaction from metal to vacancy
    type[i] = 3;
    update_list[0] = i;
    n_update_list = 1;
  }
  // update neighbor lists
  // if type_list have not all updated, keep looping
  while (n_update_list != 0){
    n_update_list = update_neighbor_check(n_update_list); // assign new list to old list
  }
}
*/
/* ----------------------------------------------------------------------
  update region_loop
  !!currently check every neighbor if it is salt region(type3)
  !!may check the effeciency in the future
  comment but keep it
------------------------------------------------------------------------- */
/*
int Appcoros::update_neighbor_check(int l)
{
// i2 = element[i]; i1 = type[i]
// i1 = 1:bulk region; =2 interface region; =3 salt region
  int num_list = 0;
  int k,m,m2,kd, kd2;
  int numneighmetal = 0;
  int n;
  int temp_update_list[200];
//check the neighbor of each lists
  for (k = 0; k<n_update_list; k++){ // loop to n_update_list
    for (m = 0; m < numneigh[update_list[k]] ; m++){ // loop of first neighbor
      kd = neighbor[update_list[k]][m];
      numneighmetal = 0; // reset to zero every loop
      //if(element[kd] == 2 && type[kd] != 3){ //set vacancy as salt region
      for (m2=0; m2<numneigh[kd]; m2++){ // loop of neighbors of neighbor
        kd2 = neighbor[kd][m2]; //site kd2 = neighbor of site kd
        if (type[kd2] == 3){ //count neighbor of kd if kd2 is salt region, then ++;
            numneighmetal++;
          }
        }
      if(element[kd] == 2 && type[kd] !=1 &&numneighmetal == 0){//if no salt region around vac, that vac become bulk region
        type[kd] =1;
        num_list++;
        //update_list[k] = kd; //store neighborlist for next check cycle
        temp_update_list[num_list-1] = kd; // temporary store neighborlist for next check cycle
      }
      if(element[kd] == 2 && type[kd] !=3 &&numneighmetal > 0){//if any salt region around vac, that vac become salt region
        type[kd] =3;
        num_list++;
        //update_list[k] = kd; //store neighborlist for next check cycle
        temp_update_list[num_list-1] = kd; // temporary store
      }
      if(element[kd] != 2 && type[kd] !=1 &&numneighmetal == 0){//if not salt region around metal, that metal become bulk region
        type[kd] =1;
      }
      if(element[kd] != 2 && type[kd] !=2 &&numneighmetal > 0){//if any salt region around metal, that metal become interface region
        type[kd] =2;
      }
  }
}
for (k = 0; k<num_list; k++){ // assign the temp store value to update_list array
  update_list[k] = temp_update_list[k];
}
return num_list;
}
*/
/* ----------------------------------------------------------------------
  update region check for surface diffusion
  comment but keep
------------------------------------------------------------------------- */
/*
int Appcoros::update_surface_diff(int i){
int num_list = 0;
int numneighmetal = 0;
int k,m,m2,kd, kd2;
int n;
int temp_update_list[200];
for (k = 0; k<n_update_list; k++){
  for (m = 0; m < numneigh[update_list[k]] ; m++){
    kd = neighbor[update_list[k]][m];
    numneighmetal = 0;
    //if(element[kd] == 2 && type[kd] != 3){ //set vacancy as salt region
    for (m2=0; m2<numneigh[kd]; m2++){
      kd2 = neighbor[kd][m2]; //site kd2 = neighbor of site kd
      if (type[kd2] == 3){ //count neighbor of kd if kd2 is salt region, then ++;
          numneighmetal++;
        }
      }
    if(element[kd] == 2 && type[kd] !=1 &&numneighmetal == 0){//if no salt region around vac, that vac become bulk region
      type[kd] =1;
      num_list++;
      //update_list[k] = kd; //store neighborlist for next check cycle
      temp_update_list[num_list-1] = kd; // temporary store
    }
    if(element[kd] == 2 && type[kd] !=3 &&numneighmetal > 0){//if any salt region around vac, that vac become salt region
      type[kd] =3;
      num_list++;
      //update_list[k] = kd; //store neighborlist for next check cycle
      temp_update_list[num_list-1] = kd; // temporary store
    }
    if(element[kd] != 2 && type[kd] !=1 &&numneighmetal == 0){//if not salt region around metal, that metal become bulk region
      type[kd] =1;
    }
    if(element[kd] != 2 && type[kd] !=2 &&numneighmetal > 0){//if any salt region around metal, that metal become interface region
      type[kd] =2;
    }
}
}
for (k = 0; k<num_list; k++){
  update_list[k] = temp_update_list[k];
}
return num_list;
}
*/
/* ----------------------------------------------------------------------
  count number of salt in lattice if i3 = 1, output to diag_coros.cpp file
------------------------------------------------------------------------- */
int Appcoros::count_salt(){
int num_salt = 0;
int i;
for (i = 0; i < nlocal; i++){
  if(potential[i] == 1){
    num_salt ++;
  }
}
return num_salt;
}
/* ----------------------------------------------------------------------
  function for salt/potential diffusion
  may be improved in count_salt part
------------------------------------------------------------------------- */
void Appcoros::potential_diff(int nsalt){
int i, j,k, kd;
int iid, jid;
int check_label;
int iwhich;
int sum = 0;
int sum_kd = 0;
//int nsalt;
double rand_i;
//nsalt = count_salt(); // count salt in lattice
  if (nsalt == 0){return;} // exception, if no salt, no diffusion

// step 1 find one randon site i which potential = 1 from all lattice(nlocal)
  iwhich = -1;
  check_label= 0;
  while(iwhich < 0 || iwhich >= nlocal || check_label !=1){
    rand_i = rancoros->uniform(); // create random float value from 0~1
    iwhich = static_cast<int> (nlocal*rand_i);
    if (potential[iwhich] == 1){
      check_label = 1;
      iid = iwhich;
    }
  }
// step 2 find a randon neighbor of site iid
// create neighbor_list of site iid
//int neighbor_list[numneigh[iid]]; // length of neighbor
int neighbor_list[numneigh[iid]] = { -1 };
  for (j = 0; j < numneigh[iid] ; j++){
    kd = neighbor[iid][j];
    if (potential[kd] == 0 && element[kd] == 2){
      neighbor_list[j] = kd;
      sum_kd ++; // number of effective neighbor
    }
    else{neighbor_list[j] = -1;}
    sum += neighbor_list[j];
  }
  if (sum == -1*numneigh[iid]){return;}// exception if all neighbor are salt or metal

 // another random draw from arrays, keep it tentatively

iwhich = -1;
check_label= 0;
int array_size = sizeof(neighbor_list)/sizeof(neighbor_list[0]); // length of array
if (array_size == 0){return;} //exception if array length = 0
 while ( iwhich < 0 || iwhich >=  array_size || check_label ==0 ){
   iwhich = rand() %  array_size;
   if (neighbor_list[iwhich] != -1 ){
     check_label = 1;
     jid = neighbor_list[iwhich];
   }
 }

// step 3 swap the potential of iid and jid
k = potential[iid];
potential[iid] = potential[jid];
potential[jid] = k;

// constant check, i1 =3 : infinite salt source
if (type[iid] ==3){potential[iid]=1;}
if (type[jid] ==3){potential[jid]=1;}

// need to update propensity after salt diffusion
// because I change the property of lattice site
//!! we can move the update_propensity after all potential_diff finished
//!! update whole lattice, be better is nmix is larger than "3000"
update_propensity(iid);
update_propensity(jid);
num_saltdiffusion ++; //count number of salt diffusion motion
}

/* ----------------------------------------------------------------------
salt_potential remove after reaction by LC
------------------------------------------------------------------------- */
void Appcoros::salt_remove(int i){
int k, kd, jid;
int iwhich, check_label;
int sum = 0;
  // create neighbor_list of site i which potential is 1
  int neighbor_list[numneigh[i]] = { -1 };
  for (k = 0; k < numneigh[i] ; k++){
    kd = neighbor[i][k];
    if (potential[kd] == 1){ //kd has impurity
      neighbor_list[k] = kd;
    }
    else {
      neighbor_list[k] = -1;
      sum ++;
    }
    //sum += neighbor_list[k];
  }
  if (numneigh[i] == 0){return;} // exception: no neighbor
  if (sum == numneigh[i]){return;}// exception: no impurity at neighbor

  // find random neighbor which site is salt
  int array_size = sizeof(neighbor_list)/sizeof(neighbor_list[0]);
  if (array_size == 0){return;} //exception
  iwhich = -1;
  check_label = 0;
  while (iwhich < 0 || iwhich >  array_size || check_label == 0){
    iwhich = rand() %  array_size;
    if (neighbor_list[iwhich] != -1 ){
      check_label = 1;
      // assign array id to lattice id
      jid = neighbor_list[iwhich];
    }
  }
  // check point
  //fprintf(screen,"[iwhich]: %d, neighbor_list[iwhich]: %d \n", iwhich, neighbor_list[iwhich]);

  // assign array id to lattice id
  //jid = neighbor_list[rand_index];
  // remove salt by set potential 1 -->> 0
  potential[jid] = 0;
  potential[i] = 0;

  update_propensity(i);
  update_propensity(jid);

  // LC sotroing active site i and j
  cur_i = i;
  cur_j = jid;

  nreact++; // count total number of reaction by LC
}

/* ----------------------------------------------------------------------
   check if perform salt diffusion
   input t is actual total time from app_lattice.cpp
------------------------------------------------------------------------- */
void Appcoros::check_saltdiffusion(double t)
{

  bigint nmix = 0;
  int nsalt = 0;
  double time_interval = 1e7;
  if(sectorflag == 0){
    if (t / time_interval >= temp_salt_time){
      temp_salt_time ++;

    for(int i = 0; i < nsaltdiffusion; i ++) {
      if(t/salt_bfreq[i] < 1){  fprintf(screen,"warning:t=%f t/salt_bfreq<1 \n", t);}
      if(t/salt_bfreq[i] < 1){break;}

       salt_time_new[i] = static_cast<bigint>(t/salt_bfreq[i]); // salt_bfreq = 100 as default
       nmix = salt_time_new[i] - salt_time_old[i];

       nsalt = count_salt(); // count salt in lattice

       while (nmix) {  //perform mixing nmix times
         nmix --;
         potential_diff(nsalt);
         if(nmix == 0) salt_time_old[i] = salt_time_new[i];  //update time
       }
      }
    }
  }


  if(sectorflag==1){ // for sector setting
  for(int i = 0; i < nsaltdiffusion; i ++) {
    if(t/salt_bfreq[i] < 1){  fprintf(screen,"warning:t=%f t/salt_bfreq<1 \n", t);}
    if(t/salt_bfreq[i] < 1){break;}
     salt_time_new[i] = static_cast<bigint>(t/salt_bfreq[i]); // salt_bfreq = 100 as default
     nmix = salt_time_new[i] - salt_time_old[i];

     nsalt = count_salt(); // count salt in lattice

     while (nmix) {  //perform mixing nmix times
       nmix --;
       potential_diff(nsalt);
       if(nmix == 0) salt_time_old[i] = salt_time_new[i];  //update time
    }
  }
  }

}

/* ----------------------------------------------------------------------
  grow memory for salt diffusion
------------------------------------------------------------------------- */
void Appcoros::grow_saltdiffusion()
{
  int n = nsaltdiffusion + 1;
  int m = 3;
  memory->grow(salt_bfreq,n,"app/coros:salt_bfreq");
  memory->grow(salt_time_old,n,"app/coros:salt_time_old");
  memory->grow(salt_time_new,n,"app/coros:salt_time_new");
}

/* ----------------------------------------------------------------------
  print every update propensity and barrier no matter, in "site_propensity" function
  it will only focus on the updated one.
  preventing double counting.
  format:
  rstyle: 1 propensity: 1.138301e-05 frequency: 4.470000 bond_ratio: 1.000000 barrier: 0.525000
------------------------------------------------------------------------- */
void Appcoros::barrier_print(int r, double propensity, double frequency,double ratio, double barrier){

  fprintf(screen,"rstyle: %d propensity: %e frequency: %f surface_effect: %f barrier: %f \n",r, propensity,frequency,ratio, barrier);
}

/* ----------------------------------------------------------------------
  Positive check surrounding neighbor
  and print out its surrounding neighbor atoms
------------------------------------------------------------------------- */
void Appcoros::np_check(int i, int jid){
// i is vacancy site
// jid is Ni or Cr site for be diffusion

int i_Ni=0, i_vac=0, i_Cr=0, j_Ni=0, j_vac=0, j_Cr=0;
int ki, kj;
for (int n=0; n< numneigh[i]; n++){
  ki = neighbor[i][n];
  if(element[ki] == 1 ){i_Ni ++;}
  if(element[ki] == 2 ){i_vac ++;}
  if(element[ki] == 3 ){i_Cr ++;}
}

for (int m=0; m< numneigh[jid]; m++){
  kj = neighbor[jid][m];
  if(element[kj] == 1 ){j_Ni ++;}
  if(element[kj] == 2 ){j_vac ++;}
  if(element[kj] == 3 ){j_Cr ++;}
}

  double total_eng = 0.0;
  double intrinsic_barrier = 0.0;
  double eng0i, eng0j, eng1i, eng1j; //energy before and after jump
  int estyle = engstyle;
  int iele = element[i];
  int jele = element[jid];
  intrinsic_barrier = mbarrier[element[jid]];
  eng0i = sites_energy(i,estyle); //broken bond with i initially,
  eng0j = sites_energy(jid,estyle); //broken bond with j initially
  // switch the element and recalculate the site energy
  element[i] = jele;
  element[jid] = iele;
  eng1i = sites_energy(i,estyle); //broken bond with i initially,
  eng1j = sites_energy(jid,estyle); //broken bond with j initially
  // switch back
  element[jid] = jele;
  element[i] = iele;
  //barrier = migbarrier * surface_effect + (eng_after - eng_before)/2.0;

  total_eng = (eng1i + eng1j - eng0i -eng0j)/2;
  fprintf(screen,"type of [i]: %d type of [jid]: %d \n",element[i], element[jid]);
  fprintf(screen,"i_Ni: %d i_vac: %d i_Cr: %d j_Ni: %d j_vac: %d j_Cr: %d \n", i_Ni, i_vac, i_Cr, j_Ni, j_vac, j_Cr);
  fprintf(screen,"intrinsic_barrier: %f total_eng_change: %f id_i: %d id_j: %d \n",intrinsic_barrier, total_eng, i, jid);
}

/* ----------------------------------------------------------------------
  calculating total metal energy which the particle is not VAC
  // not realistic for current model
------------------------------------------------------------------------- */
double Appcoros::total_metal_energy( )
{
  double penergy = 0.0;
  for(int i = 0; i < nlocal; i++){
  if(element[i] != VACANCY){
    penergy += sites_energy(i,engstyle);
    }
  }
/*
  if(elastic_flag) {
    for(int j = 0; j < nlocal; j++) {
      int jtype = element[j];
      penergy += elastic_energy(j,jtype);
    }
  }
*/
  return penergy/2.0;
}

/* ----------------------------------------------------------------------
  stop criteria based on number of corrosion event
  first set 50%
------------------------------------------------------------------------- */
int Appcoros::KMC_stop(){
  if (KMC_stop_flag==0) return 0;
  double Cr_num_threshold = threshold_Cr * 0.01 * total_Cr;

  double j = (float) nreact;
  if (j >= Cr_num_threshold){
    //fprintf(screen,"stop by reaching nreact\n");
    return 1;
  }
  return 0;
}
/* ----------------------------------------------------------------------
vac monomer count for all lattice --> call by diag_coros.cpp
potential[i] =0 and potential[jd] = 0 then sum_monomer++;
this function count by i3 value not i2 and only works for pure vac
------------------------------------------------------------------------- */
int Appcoros::vac_monomer_count(){

  int sum_monomer = 0;
  int monomer_flag;
  int jd;
  for(int i = 0; i < nlocal; i++){  // loop all lattice
    monomer_flag = 0;
    if (potential[i] == 0){
      for (int k = 0; k < numneigh[i];k++){
          jd = neighbor[i][k];
          if (potential[jd] == 0){
            monomer_flag = 0;
            break;
          }else{monomer_flag = 1;}
      }
      if (monomer_flag == 1){sum_monomer++;}
    }
  }
return sum_monomer;
}

/* ----------------------------------------------------------------------
all monomer count for all lattice --> call by diag_coros.cpp
this function count by i2 value not i3
------------------------------------------------------------------------- */
void Appcoros::monomer_count(){

  int sum_monomer = 0;
  int monomer_flag;
  int jd;

  // reset monomers array
  for(int i = 1; i < nelement; i++) {
    monomers[i] = 0;
  }
  for(int element_i = 1; element_i < nelement; element_i++) { // loop through all element
    for(int lattice_i = 0; lattice_i < nlocal; lattice_i++) {  // loop all lattice
      monomer_flag = 0;
      if (element[lattice_i] == element_i) {
        for (int k = 0; k < numneigh[lattice_i];k++){
          jd = neighbor[lattice_i][k];
            if (element[jd] == element_i){
              monomer_flag = 0;
              break;
              }else{monomer_flag = 1;}
          }
          if (monomer_flag == 1){monomers[element_i]++;}
      }
    }
}

}

/* ----------------------------------------------------------------------
  The concentration_field func and time_averaged_concentration func is extract from
  app_seg.cpp which is originated developed by Yongfeng
  Integrate c*t at each site for fractional occupancy over time
  the count for time_average_concentration and site_concentration is in this function.
  but separate into different func due to different d_time when average compute

  // update dt_site with better efficiency
  // the idea of site concentration are the same
  // but change the implimentation, update the active site only
  // other site will be update at end of sector
------------------------------------------------------------------------- */
void Appcoros::concentration_field(double dt)
{
  //dt_new += dt; // update time interval for time_average_concentration

  //dt_i3_site_new += dt;

  // for(int i = 0; i < nlocal; i++) {
  //    ct_new[element[i]] += dt; // total concentration
  //    ct_site_new[element[i]][i] += dt; // site concentration per kmc step
  //    i3_site_new[potential[i]][i] += dt;
  // }

  // total occupancy on whole lattice, different (shorter) time interval from site concentration
  ct_time += dt;
  ct_new[1] += total_Ni * dt;
  ct_new[2] += total_vac * dt;
  ct_new[3] += total_Cr * dt;

  //active_vac_new share same time interval with ct_time

  // LC, new site concentration calculation, testing
  //i and j are the active site that do site event
  // this works for diffusion and reaction
  dt_site_c_new += dt;       // total time interval in 1 sector
  //site_time_interval += dt;  // total time interval across sectors, total time for site concentration
  time_sector += dt;         // current time point in sector;


  double t_interval = 0.0;
  t_interval = (time_sector - dt ) - ct_site_temp_t[cur_i];
  ct_site_new[ct_site_temp_i2[cur_i]][cur_i] += t_interval ;  // occupancy
  site_time[cur_i] += t_interval;
  // LC test
  if(t_interval<0){
    fprintf(screen, "t_interval:%f, time_sector:%f, dt:%f, ct_site_temp_t[cur_i]:%f\n", t_interval, time_sector, dt, ct_site_temp_t[cur_i]);
  }


  if (ct_site_temp_i2[cur_i]==2) {active_vac_new += t_interval;} // active vac
  //active_vac_new += t_interval;

  t_interval = (time_sector - dt ) - ct_site_temp_t[cur_j];
  ct_site_new[ct_site_temp_i2[cur_j]][cur_j] += t_interval ;
  site_time[cur_j] += t_interval;

  // LC test
  if(t_interval<0){
    fprintf(screen, "t_interval:%f, time_sector:%f, dt:%f, ct_site_temp_t[cur_i]:%f\n", t_interval, time_sector, dt, ct_site_temp_t[cur_j]);
  }

  if (ct_site_temp_i2[cur_j]==2) {active_vac_new += t_interval;} // active vac

  //fprintf(screen, "time_sector:%f, dt:%f\n", time_sector, dt);
  //update temp time, time format
  ct_site_temp_t[cur_i] = (time_sector - dt );  // update
  ct_site_temp_t[cur_j] = (time_sector - dt );
  // LC test
  if(time_sector-dt<0){
    fprintf(screen, "time_sector-dt:%f\n", (time_sector-dt));
  }
  // update, element/i2 type
  ct_site_temp_i2[cur_i] = element[cur_i]; // update
  ct_site_temp_i2[cur_j] = element[cur_j];

  return;
}

/* ----------------------------------------------------------------------
concentration_field for global algorithm
------------------------------------------------------------------------- */
void Appcoros::concentration_field_global(double dt)
{
  // total occupancy on whole lattice, different (shorter) time interval from site concentration
  ct_time += dt;
  ct_new[1] += total_Ni * dt;
  ct_new[2] += total_vac * dt;
  ct_new[3] += total_Cr * dt;
  dt_site_c_new += dt;
  time_sector += dt;

  // parameter for diag version
  time_sector_diag += dt;
  dt_site_c_new_diag += dt;

  double t_interval = 0.0;
  double t_interval_diag = 0.0;
  t_interval = (time_sector - dt ) - ct_site_temp_t[cur_i];
  ct_site_new[ct_site_temp_i2[cur_i]][cur_i] += t_interval ;  // occupancy
  site_time[cur_i] += t_interval;
  // diag version:
  t_interval_diag = (time_sector_diag - dt ) - ct_site_temp_t_diag[cur_i];
  ct_site_new_diag[ct_site_temp_i2_diag[cur_i]][cur_i] += t_interval_diag ;  // occupancy
  site_time_diag[cur_i] += t_interval_diag;
  if (ct_site_temp_i2[cur_i]==2) {active_vac_new += t_interval;} // active vac


  t_interval = (time_sector - dt ) - ct_site_temp_t[cur_j];
  ct_site_new[ct_site_temp_i2[cur_j]][cur_j] += t_interval ;
  site_time[cur_j] += t_interval;
  // diag version:
  t_interval_diag = (time_sector_diag - dt ) - ct_site_temp_t_diag[cur_j];
  ct_site_new_diag[ct_site_temp_i2_diag[cur_j]][cur_j] += t_interval_diag ;  // occupancy
  site_time_diag[cur_j] += t_interval_diag;

  if (ct_site_temp_i2[cur_j]==2) {active_vac_new += t_interval;} // active vac

  //update temp time, time format
  ct_site_temp_t[cur_i] = (time_sector - dt );
  ct_site_temp_t[cur_j] = (time_sector - dt );

  // diag version update
  ct_site_temp_t_diag[cur_i] = (time_sector_diag - dt );
  ct_site_temp_t_diag[cur_j] = (time_sector_diag - dt );

  // update, element/i2 type
  ct_site_temp_i2[cur_i] = element[cur_i]; // update
  ct_site_temp_i2[cur_j] = element[cur_j];

  // diag version update
  ct_site_temp_i2_diag[cur_i] = element[cur_i]; // update
  ct_site_temp_i2_diag[cur_j] = element[cur_j];

}
/* ----------------------------------------------------------------------
update site concentration
for sector setting
by LC
------------------------------------------------------------------------- */
void Appcoros::site_concentration_calc(int iset, int nset){
  // iset: nth sector
  // for n = 4, iset = 0, 1, 2, 3
  int j;
  double t_interval = 0.0;
  int first = 0 + iset * nlocal/nset;
  int last = (iset+1) * nlocal/nset;
  for (j=first;j<last;j++){
    t_interval = dt_site_c_new - ct_site_temp_t[j];

    ct_site_new[ct_site_temp_i2[j]][j] += t_interval ;
    // site experienced time; !! by LC
    site_time[j] += t_interval;
    //ct_site_temp_t[j] = 0.0; // reset
  }

  fill_n(ct_site_temp_t, nlocal, 0.0); // set whole array to 0
  site_time_interval[iset] += time_sector;  // add site_time_interval for each sector

    // for (j=0;j<nlocal;j++){
    //   t_interval = dt_site_c_new - ct_site_temp_t[j];
    //
    //   // LC test
    //   if(t_interval<0){
    //     fprintf(screen, "t_interval:%f, dt_site_c_new:%f, ct_site_temp_t[j]:%f\n", t_interval, dt_site_c_new, ct_site_temp_t[j]);
    //   }
    //
    //   ct_site_new[ct_site_temp_i2[j]][j] += t_interval ;
    //   ct_site_temp_t[j] = 0.0; // reset
    //
    // }
    dt_site_c_new = 0.0; //reset
    time_sector = 0.0;  // reset
}

/* ----------------------------------------------------------------------
  update time averaged total concentration concentrations only
  this function runs after every sector loop
------------------------------------------------------------------------- */
void Appcoros::time_averaged_concentration()
{
  // if(dt_new <= 0){
  //   return;
  // }
  // for(int i = 1; i <= nelement; i++) { // ct = c*t / dt
  //    ct[i] = ct_new[i]/dt_new/nlocal; //time and spatial average
  //    ct_new[i] = 0.0; //start recounting
  //  }
     // !!site_concentration part was separated to ct_site_extract!!
     // for(int j = 0; j < nlocal; j++) {
	   //    //disp[ndiff+i][j] = ct_site[i][j]/dt_new; //time average
     //    ct_site[i][j] = ct_site_new[i][j]/dt_new;
	   //    ct_site_new[i][j] = 0.0; //start recounting
     // }

// LC
  // int i,j;
  //  for(i = 1; i <= nelement; i++) {
  //     for(j = 0; j < nlocal; j++) {
  //        ct_site[i][j] = ct_site_new[i][j]/site_time_interval;
  //        ct_site_new[i][j] = 0.0; //recounting
  //     }
  //  }
     // LC
     // for (int j = 0; j<nlocal; j++){
     //   ct_site_temp_t[j] = 0.0; // reset
     // }
     //dt_site_c_new = 0.0;

  //dt_new = 0.0;
}

/* ----------------------------------------------------------------------
ct_site_extract
compute site_concentration and average
return ct_site **array called by app_lattice <-- dump_text
use ctNi, ctVac, ctCu as input dump
------------------------------------------------------------------------- */
double **Appcoros::ct_site_extract(){

  // for global setting
  if (sectorflag ==0){
    if(dt_site_c_new <= 0){  // for first dump
      return ct_site;
    }
    int i ,j;
    double t_interval;
    for (j=0;j<nlocal;j++){
      t_interval = dt_site_c_new - ct_site_temp_t[j];
      ct_site_new[ct_site_temp_i2[j]][j] += t_interval ;
      // site experienced time; !! by LC
      site_time[j] += t_interval;
  }
    fill_n(ct_site_temp_t, nlocal, 0.0); // set whole array to 0

    dt_site_c_new = 0.0; //reset
    time_sector = 0.0;  // reset

    for(i = 1; i <= nelement; i++) {
       for(j = 0; j < nlocal; j++) {
          ct_site[i][j] = ct_site_new[i][j]/site_time[j];
          ct_site_new[i][j] = 0.0; //recounting
       }
    }
      for (j=0; j<nlocal; j++){
        site_time[j] = 0.0;
      }
    }
  return ct_site;

  // for sector setting
  if (sectorflag ==1){
  int i,j;
  if(site_time_interval[0] <= 0){
    return ct_site;
  }


   for(i = 1; i <= nelement; i++) {
      for(j = 0; j < nlocal; j++) {
         ct_site[i][j] = ct_site_new[i][j]/site_time[j];
         ct_site_new[i][j] = 0.0; //recounting
      }
   }

   for (j=0; j<8;j++){
   site_time_interval[j] = 0.0;
 }

 for (j=0; j<nlocal; j++){
   site_time[j] = 0.0;
 }
}
  // LC summation rest of them when compute


  // for(int i = 1; i <= nelement; i++) {
  //    for(int j = 0; j < nlocal; j++) {
  //       ct_site[i][j] = ct_site_new[i][j]/dt_site_c_new;
	//       ct_site_new[i][j] = 0.0; //recounting
  //    }
  // }
  // dt_site_c_new = 0.0;
  // // LC
  // for (int j = 0; j<nlocal; j++){
  //   ct_site_temp_t[j] = 0.0; // reset
  // }



//// special case analysis
  // neighbor check and print if ct_site[vac][j] > 0 but its neighbor ct_site[vac][neighbor] is 0
  // special case that if vac exist, it should diffuse to neighbor to create ct_site
  // shoud comment for normal use
  // for(int j = 0; j < nlocal; j++) {
  //   check_flag =0;
  //   if (ct_site[2][j]>0){
  //   for (int k = 0; k < numneigh[j];k++){
  //     jd = neighbor[j][k];
  //     //i ==2 == vac
  //     if (ct_site[2][jd] >0.0){
  //     check_flag = check_flag + 1;}
  //     else{check_flag = check_flag + 0;}
  //     }
  //   }
  //   if (check_flag ==0&&ct_site[2][j]>0) {fprintf(screen,"check point\n");
  //     } // if check point
  //   }
//// special case analysis

  return ct_site;
}

/* ----------------------------------------------------------------------
ct_extract
compute site_concentration and average
return ct_site **array called by app_lattice <-- dump_text
use ctNi, ctVac, ctCu as input dump
------------------------------------------------------------------------- */
double *Appcoros::ct_extract(){
  int i,j;
if(ct_time <= 0){
  return ct;
}

  for (i=1; i<= nelement; i++){
    ct[i] = ct_new[i]/ct_time/nlocal;
    ct_new[i] = 0.0; //reset
  }
ct_time = 0.0; // reset

return ct;
}
/* ----------------------------------------------------------------------
active_vac_extract
from diag_
------------------------------------------------------------------------- */
double Appcoros::active_vac_extract(){
  int i,j;
if(ct_time <= 0){
  return active_vac;
}

active_vac = active_vac_new / ct_time/nlocal;
active_vac_new = 0.0;

return active_vac;
}
/* ----------------------------------------------------------------------
site concentration of i3 ==0,1
output in dump file per site
------------------------------------------------------------------------- */
double **Appcoros::i3_site_extract(){
//output: number of pure_vac as integer in stats/diag
//int pure_vac;

if(dt_i3_site_new <= 0){
  return i3_site;
}
// if(concentrationflag==0){ //excpetion for no concentrationflag
//   return i3_site;
// }
// count site concentration of i3==0 and i3==1
for(int i = 0; i < 2; i++) {
  for(int j = 0; j < nlocal; j++){
    i3_site[i][j] = i3_site_new[i][j]/dt_i3_site_new;
    i3_site_new[i][j] = 0.0;
  }
}
dt_i3_site_new = 0;

return i3_site;
}
/* ----------------------------------------------------------------------
output temp into in dump file per site
------------------------------------------------------------------------- */
double *Appcoros::temp_t_extract(){

return ct_site_temp_t;
}
/* ----------------------------------------------------------------------
output temp into in dump file per site
------------------------------------------------------------------------- */
int *Appcoros::temp_i2_extract(){

return ct_site_temp_i2;
}
/* ----------------------------------------------------------------------
metal_pure_vac approximation in lattice
this function is count the pure vac in metal region in lattice without knowing the boundary
Concept:
condition 1:count number of i3 ==0 and its site_c of (i3==1) ==0
condition 2: neighbor of i must not be impurity(i3=1)
for each time interval,this particle must be i3==0 at first and site_c of impurity ==0,
if site_concentration of (i3 ==1) >0, which means this position is pollute by impurity diffusion

------------------------------------------------------------------------- */
int Appcoros::metal_pure_vac_approxi()
{
  int metal_pure_vac =0;
  int jd;
  int sum_neigh;
  //if potential[i] == 0 and ct_site[1][i] ==0;

  //if no concentration flag no cite concentration, the value is wrong
  // if (concentrationflag==0){
  //   return 0;
  // }

  for(int i =0; i<nlocal;i++){
    if (potential[i] == 0 && i3_site[1][i] ==0.0){
      sum_neigh=0;
      for (int k = 0; k < numneigh[i];k++){
        jd = neighbor[i][k];
        if (potential[jd]==1){
          sum_neigh++;
        }
      }
      if (sum_neigh==0){
        metal_pure_vac++;
        //fprintf(screen,"metal_pure_vac_id:%d\n", i); // print out metal_pure_vac id for check
      }
    }
  }

  return metal_pure_vac;
}
/* ----------------------------------------------------------------------
dump all event_list per dump time, output as txt file by LC
this function is called from app.lattice, can be improve in the future
dump for every sector loop
------------------------------------------------------------------------- */
void Appcoros::dump_event(double dt_kmc)
{
  //if(dump_event_flag != 1) return;

  char temp[10];
  string dumpfilename = ".dump_eventlist";
  sprintf(temp, "%d", dump_index);
  dumpfilename = temp + dumpfilename;
  //create/open file
  ofstream file(dumpfilename);
  string header = "index event_id style which jpartner element[j] next propensity\n";
  file << "number of event: " << nevents << "\n";
  file << "dt_sector: " << dt_kmc << "\n";
  file << header;
  //print every event list

  // dump based on event id, may include unwanted event
  // for (int i = 0; i < nevents; i++){
  //   if(element[events[i].jpartner]!=2){
  //   file << i << " "
  //        << events[i].style << " "
  //        << events[i].which << " "
  //        << events[i].jpartner << " "
  //        << element[events[i].jpartner] << " "
  //        << events[i].next << " "
  //        << events[i].propensity << " \n";
  //      }
  // }


int event_id;
int idx=0;
//dump based on atom site
  for (int i = 0; i < nlocal; i++){
    event_id = firstevent[i]; // first event id of own sites i
    while(event_id!=-1){
    file << idx << " "
         << event_id << " "
         << events[event_id].style << " "
         << events[event_id].which << " "
         << events[event_id].jpartner << " "
         << element[events[event_id].jpartner] << " "
         << events[event_id].next << " "
         << events[event_id].propensity << " \n";
    idx++;
    event_id = events[event_id].next;
       }
  }

  dump_index++;
  file.close();
}

/* ----------------------------------------------------------------------
dump all event_list in one file
this function dump the event information during add_event
------------------------------------------------------------------------- */
void Appcoros::dump_event_all(int style, int which, int jpartner, int next, double propensity)
{
string filename = "all.dump_eventlist";
if(dump_event_all_file == 0){
  ofstream file(filename);
  string header = "index style which jpartner next propensity\n";
  file << header;
  dump_event_all_file = 1;
  file.close();
}

//freopen(file, filename);
ofstream file (filename, ios_base::app);
file << dump_event_all_index << " "
     << style << " "
     << which << " "
     << jpartner << " "
     << next << " "
     << propensity << " \n";
dump_event_all_index ++;
file.close();
}

/* ----------------------------------------------------------------------
MSD calculation
aid is atomic id which is change with event
------------------------------------------------------------------------- */
void Appcoros::MSD_calculation()
{

  //fprintf(screen, "x,y,z: %f, %f, %f\n", xyz[i][0],xyz[i][1],xyz[i][2]);

//double Ni_disp;
//double vac_disp;
double x_disp, y_disp, z_disp;
//double sd = 0.0; //square displacement
int idx, current_id, initial_id;
for(int i = 1; i < nelement; i++) {
  sd[i] = 0.0; // declare at header
}
//idx: position id
//aid: current atom id
//element[current_id]  //i2 of atom
for ( idx = 0 ;idx < nlocal; idx++){
  initial_id = idx;
  current_id = aid[idx];
  x_disp = fabs(xyz[current_id][0] - xyz[initial_id][0]);
  y_disp = fabs(xyz[current_id][1] - xyz[initial_id][1]);
  z_disp = fabs(xyz[current_id][2] - xyz[initial_id][2]);
  sd[element[current_id]] = sd[element[current_id]] + x_disp*x_disp + y_disp*y_disp + z_disp*z_disp;
  }

  //fprintf(screen,"sd[1]:%f\n", sd[1]);
  //fprintf(screen,"sd[2]:%f\n", sd[2]);
  // if (element[current_id]==1){//Ni msd
  //   sd[i] = sd[i] + x_disp*x_disp + y_disp*y_disp + z_disp*z_disp;  // summation of square displacement
  //
  //   }

}

/* ----------------------------------------------------------------------
coros info output
This function will be called when corrosion happend
this will record the time, Cr position, corrosion rate
------------------------------------------------------------------------- */
void Appcoros::corrosion_info_output(int i)
{
// i in Cr position id
// dt as the coros time for corrosion
  //fprintf(screen, "x,y,z: %f, %f, %f\n", xyz[i][0],xyz[i][1],xyz[i][2]);
  //fprintf(screen, "time1: %f \n", time);  // --> time point for corrosion event
  //fprintf(screen, "dt_step1: %f \n", dt_step);  // --> dt_step event for corrosion
  coros_time[coros_index] = time;
  coros_id[coros_index] = i;  // position id ~= position

  coros_index ++;
}

/* ----------------------------------------------------------------------
coros info dump:
transform id to position and dump array to one log file
------------------------------------------------------------------------- */
void Appcoros::corrosion_info_dump()
{

  char temp[10];
  string dumpfilename = "coros_info";
  //fprintf(screen, "check \n" );
  //create/open file
  ofstream file(dumpfilename);
  string header = "index time x y z\n";
  file << header;
  //dump based on atom site
    for (int i = 0; i < nlocal; i++){
      //fprintf(screen, "coros_id%d\n",coros_id[i] );
      if(coros_id[i]==-1) break;
      file << i << " "
           << coros_time[i] << " "
           << xyz[coros_id[i]][0] << " "  // x
           << xyz[coros_id[i]][1] << " "  // y
           << xyz[coros_id[i]][2] << " \n";  // z
    }

    file.close();
}

/* ----------------------------------------------------------------------
diag_site concentration
this site concentration is compute to bulk_cvac
this only work for global(no sector)
it may affect the efficiency -> need to check
------------------------------------------------------------------------- */
void Appcoros::site_concentration_diag()
{
  // for global setting

    int i ,j;
    double t_interval_diag;
    for (j=0;j<nlocal;j++){
      t_interval_diag = dt_site_c_new_diag - ct_site_temp_t_diag[j];
      ct_site_new_diag[ct_site_temp_i2_diag[j]][j] += t_interval_diag ;
      site_time_diag[j] += t_interval_diag;
  }
    fill_n(ct_site_temp_t_diag, nlocal, 0.0); // set whole array to 0

    dt_site_c_new_diag = 0.0; //reset
    time_sector_diag = 0.0;  // reset

    for(i = 1; i <= nelement; i++) {
       for(j = 0; j < nlocal; j++) {
          ct_site_diag[i][j] = ct_site_new_diag[i][j]/site_time_diag[j];
          ct_site_new_diag[i][j] = 0.0; //recounting
       }
    }
      for (j=0; j<nlocal; j++){
        site_time_diag[j] = 0.0;
      }

  //return ct_site_diag;

}

/* ----------------------------------------------------------------------
bulk_cvac
this function call extract site concentration to compute bulk cvac
by determining the interface
this process will flatten/average the site concentration into 1D for process
------------------------------------------------------------------------- */
double Appcoros::bulk_cvac()
{
  if(dt_site_c_new_diag > 0){  // for first dump
    site_concentration_diag();  // update site concentration first;
  }

  double bulk_cvac;
  int i, z;
  int temp_Cr_pos = -1;
  // reset profile
  for (i = 0; i < z_height; i++)
  {
    zprofile[i] = 0.0;
  }

  for (i = 0; i < nlocal; i++)
  {
    // z axis: xyz[i][2]
    //ct_site[i][j] : i for element, j for site_id
    //ct_site[2][i]
    z = int(xyz[i][2]*2);  // transfer z axis to 1D profile id
    zprofile[z] += ct_site_diag[2][i];
    if (element[i] == 1 && z > temp_Cr_pos){
      temp_Cr_pos = z;
    }

  }


  // find interface point from profile by mininum point for Cvac
  // create cvac profile along z
  for (i = 0; i < z_height; i++)
  {
    zprofile[i] = zprofile[i] / xy_sum ;
    //fprintf(screen, "zprofile[%d]: %e\n",i,  zprofile[i]);
  }
  // find interface

  // method 1 find mininum site vac concentration
  int min_pos=-1;
  double min = 1.0;
  for (i = int(z_height*0.8); i < z_height; i++){
    if (zprofile[i] < min && zprofile[i] >0.0){
      min = zprofile[i];
      min_pos = i;
    }
  }
  fprintf(screen, "min_pos:%d, min:%e\n", min_pos, min );

  // method 2: find maximum Cr site as interface
  temp_Cr_pos = temp_Cr_pos -1; // LC test
  fprintf(screen, "Cr_pos:%d\n", temp_Cr_pos);
  // summation all site vac concentration to the height of interface as bulk Cvac
  int temp_sum=0;
  for (i = 0; i < temp_Cr_pos; i++)
  {
    bulk_cvac =+ zprofile[i];
    temp_sum ++;
  }

  bulk_cvac = bulk_cvac / temp_sum;
  fprintf(screen, "bulk_cvac %e\n", bulk_cvac);
  return bulk_cvac;
}
