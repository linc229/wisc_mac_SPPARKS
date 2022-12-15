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
   This application does radiation segregation simulations based on Soisson 2006.
   Contributer: Yongfeng Zhang, yzhang2446@wisc.edu
------------------------------------------------------------------------- */

/* Things to be done
 * 1. Update accelerated diffusion
 * (Done) 2. Correct diffusion calculation
 * (Done) 3. Add corrections to energy and barrier for different interstitials
 * (Done) 4. Enable ct_site calculation to show the spatial distribution of concentrations
 * 5. Add full references to models
 * 6. Update reaction
*/


#include "math.h"
#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "app_seg.h"
#include "domain.h"
#include "solve.h"
#include "random_park.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"

using namespace SPPARKS_NS;

enum{NOOP,BCC,NBCC,WALL};                          // all sites are on lattice, specicial site types (e.g., sinks) can be added
enum{VAC=0,INT,CE1,CE2,CE3,CE4,CE5,CE6,CE7,CE8};   // CE: chemical element

#define DELTAEVENT 100000
#define MAX2NN 6 // max # of 2NN for both FCC and BCC lattice
#define BIGNUMBER 1e18 // define a big number

/* ---------------------------------------------------------------------- */

AppSeg::AppSeg(SPPARKS *spk, int narg, char **arg) :
  AppLattice(spk,narg,arg)
{
  ninteger = 5; // first for lattice type,second for element, 3rd for SIA type, 4th and 5th dumbbell atoms
  ndouble = 0;
  delpropensity = 2;
  delevent = 1;
  allow_kmc = 1;
  allow_rejection = 0;

  engstyle = 0; //1 for 1NN interaction, 2 for 2NN interaction; default 0
  diffusionflag = 0; //flag for MSD calculations, 1 for yes, 0 for no; default 0
  concentrationflag = 0; //flag for concentration calculation
  seg_flag = 0;
  eisink_flag = 0;

  if (narg < 1) error->all(FLERR,"Illegal app_style command");
  if (narg >= 2) engstyle = atoi(arg[1]);
  if (engstyle == 2) delpropensity += 1;// increase delpropensity for 2NN interaction
  if (narg >= 3) diffusionflag = atoi(arg[2]);
  // darray 1-4 for msd if activated, followed by concentrations if activated, needs the initial atomic id, defined as aid
  if (diffusionflag > 0) {ninteger++; ndouble += 4;}
  if (narg >= 4) concentrationflag = atoi(arg[3]);
  // calculate concentration fields for each elements, so that concentrationflag = nelement
  if (concentrationflag) {
     ndiff = 4*diffusionflag;
     ndouble += concentrationflag+1;
  }

  create_arrays();

  firsttime = 1;
  esites = NULL;
  echeck = NULL;
  events = NULL;
  maxevent = 0;
  firstevent = NULL;

  // set random number generator
  ranseg = new RandomPark(ranmaster->uniform());
  double seed = ranmaster->uniform();
  ranseg->reset(seed,me,100);

  // flags for bond interactions
  ebond1 = NULL;
  ebond2 = NULL;
  mbarrier = NULL;
  hcount = NULL; //numner of vacancy switch events
  //nn1flag = nn2flag = barrierflag = 0; //flags for bond energy and migration barriers
  nn1flag = nn2flag = barrierflag = 0; //flags for bond energy and migration barriers
  mfpflag = 0; // Peng: mfpflag

  // flags and parameters for sinks, dislocations, reactions
  sink_flag = elastic_flag = moduli_flag = dislocation_flag = reaction_flag = acceleration_flag = 0; //flags for sink dislocation and vacancy
  nsink = ndislocation = nreaction = nballistic = ntrap = 0;

  // arrays for dislocations
  dislocation_type = line_vector = nsegment = NULL;
  burgers = xdislocation = NULL;
  dislocation_radius = NULL;

  // arrays for sinks
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
  bfreq = NULL;
  time_old = time_new = NULL;
  min_bfreq = BIGNUMBER;

  // arrays for frenkel pair productionr
  fpfreq = 0;
  fp_old = fp_new = 0;
  min_fpfreq = BIGNUMBER;

  // 2NN neigbor information
  numneigh2 = NULL;
  neighbor2 = NULL;

  // dumbbell
  edumbbell = NULL;
  emdumbbell = NULL;
  vdumbbell = NULL;
  nsia = NULL;

  // ris
  ris_type = NULL;
  ris_ci = ris_total = NULL;

  // short range order
  total_neighbor = NULL;
  sro = NULL;
}

/* ---------------------------------------------------------------------- */

AppSeg::~AppSeg()
{
  delete [] esites;
  delete [] echeck;
  delete ranseg;

  memory->sfree(events);
  memory->destroy(firstevent);
  memory->destroy(ebond1);
  memory->destroy(ebond2);
  memory->destroy(mbarrier);
  memory->destroy(nsites_local);
  memory->destroy(edumbbell);
  memory->destroy(vdumbbell);
  memory->destroy(hcount);

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
  }


  if (mfpflag) {// memory use for accelerate simulation (Peng: mfpflag)
    memory->destroy(rhmfp);
    memory->destroy(nhmfp);
    memory->destroy(mfp);
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

}

/* ----------------------------------------------------------------------
  setup bond energies for 1NN, 2NN, and SPs
------------------------------------------------------------------------- */

void AppSeg::input_app(char *command, int narg, char **arg)
{
  int i,j,ibond;
  int nlattice = nlocal + nghost;

  // 1NN bond energy taken in the order of 11 12 ... 1N; 22 ... 2N; ...; NN
  if (strcmp(command,"ebond1") == 0) {

    if (narg < 3) error->all(FLERR,"Illegal ebond1 command");
    nelement = atoi(arg[0]);   // num of elements

    memory->create(nsites_local,nelement,"app/seg:nsites_local");
    memory->create(total_neighbor,nelement,"app/seg:total_neighbor"); //total number of neighbors of all type i atoms; needed for sro calculations
    memory->create(sro,nelement,nelement,"app/seg:sro"); //short range order matrix
    memory->create(ci,nelement,"app/seg:ci"); //static concentration based on current configuration
    memory->create(ct,nelement,"app/seg:ct"); //time averaged concentration based on the fractional occupation at each site
    memory->create(ct_new,nelement,"app/seg:ct_new"); //time averaged concentration
    memory->create(ebond1,nelement,nelement,"app/seg:ebond1"); // 1NN bond energy
    memory->create(hcount,nlattice,"app/erad:hcount");
    if(diffusionflag) memory->create(Lij,nelement,nelement,"app/seg:Lij"); //Onsager coefficient
    if(concentrationflag) memory->create(ct_site,nelement,nlocal,"app/seg:ct_site"); //site coefficient

    if(narg != nelement*(nelement+1)/2+1) error->all(FLERR,"Illegal ebond1 command");
    nn1flag = 1;

    for (i = 0; i < nelement; i++ ) {
      for (j = i; j < nelement; j++ ) {
        ibond = ibonde(i+1,j+1,nelement);
        ebond1[i][j] = atof(arg[ibond]);
        if(j > i) ebond1[j][i] = ebond1[i][j];
      }
    }
   }

  // 2NN bond energy taken in the order of 11 12 ... 1N; 22 ... 2N; ...; NN
  else if (strcmp(command,"ebond2") == 0) {

    nelement = atoi(arg[0]);   // num of elements, just to be consistent with ebond1
    memory->create(ebond2,nelement,nelement,"app/seg:ebond2");
    if(narg != nelement*(nelement+1)/2+1) error->all(FLERR,"Illegal ebond2 command");

    nn2flag = 1;
    for (i = 0; i < nelement; i++ ) {
      for (j = i; j < nelement; j++ ) {
        ibond = ibonde(i+1,j+1,nelement);
        ebond2[i][j] = atof(arg[ibond]);
        if (j > i) ebond2[j][i] = ebond2[i][j];
      }
    }
  }

  // formation energy of different types of SIAs defined by the two dumbbell atoms
  // Note the total element include V and I, which are "needed" here but useless
  else if (strcmp(command,"edumbbell") == 0) {

    memory->create(edumbbell,nelement,nelement,"app/seg:edumbbell");
    memory->create(emdumbbell,nelement,nelement,"app/seg:emdumbbell");
    memory->create(fdumbbell,nelement,nelement,"app/seg:fdumbbell");
    number_sia = (nelement-2)*(nelement-1)/2;
    memory->create(nsia,number_sia,"app/seg:nsia");
    for (i = 0; i < number_sia; i++ ) {nsia[i] = 0;}

    if(narg != nelement*(nelement+1)/2+1) error->all(FLERR,"Illegal edumbbell command");

    nn2flag = 1;
    for (i = 0; i < nelement; i++ ) {
      for (j = i; j < nelement; j++ ) {
        ibond = ibonde(i+1,j+1,nelement);
        edumbbell[i][j] = atof(arg[ibond]);
        if (j > i) edumbbell[j][i] = edumbbell[i][j];
      }
    }
  }

  // sia barrier correction based the "type" of sias
  else if (strcmp(command,"emdumbbell") == 0) {
    if(narg != nelement*(nelement+1)/2+1) error->all(FLERR,"Illegal emdumbbell command");

    nn2flag = 1;

    for (i = 0; i < nelement; i++ ) {
      for (j = i; j < nelement; j++ ) {
        ibond = ibonde(i+1,j+1,nelement);
        emdumbbell[i][j] = atof(arg[ibond]);
        if (j > i) emdumbbell[j][i] = emdumbbell[i][j];
      }
    }
  }

  // types of SIAs and the dumbbell vectors.
  else if (strcmp(command, "dumbbell") ==0) {
    ndumbbell = atoi(arg[0]); // number of dumbbell types
    memory->create(vdumbbell,ndumbbell,3,"app/seg:vdumbbell");
    for (i = 0; i < ndumbbell; i++) {
        vdumbbell[i][0] = atof(arg[i*3])+1;
        vdumbbell[i][1] = atof(arg[i*3])+2;
        vdumbbell[i][2] = atof(arg[i*3])+3;
    }
  }

  // migration barriers for each element. For V this is the migration barrier is set to be 0 (given by the elements).
  else if (strcmp(command, "migbarrier") ==0) {

    if (narg < nelement) error->all(FLERR,"Illegal migbarrier command, a barrier is needed for each element");
    barrierflag = 1;
    memory->create(mbarrier,nelement,"app/seg:mbarrier");

    for (i = 0; i < narg; i++) {
        mbarrier[i] = atof(arg[i]);
    }
  }

  // list of elements that interact with vacancy and interstitials.
  // default 0 meaning no trapping; 1 means trapping. Need for all elements if defined
  else if (strcmp(command, "trap") ==0) {

    if (narg < nelement) error->all(FLERR,"Illegal trap command, a trap flag is needed for each element");
    for (i = 0; i < nelement; i++) {
        trap[i] = atof(arg[i]);
    }
  }

  // elastic moduli for stress calculation, only for cubic system now
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
      memory->create(stress,nlattice,6,"app/seg:stress");
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

  // elastic interaction with stress field via relaxed volume; can be extended for dipole tensor later
  else if (strcmp(command, "elastic_interaction") ==0) {

    elastic_flag = 1;
    if(narg < 2) error->all(FLERR,"illegal elastic interaction command");

    int iarg = narg/2;
    for (i = 0; i < iarg; i++) {
      int ei = i*2;  //element
      evol[ei] = atof(arg[i*2+1]); //relaxation volume
    }
  }

  // define sinks to defects and element, one sink each line
  // both sink and sink_interaction commends need to be defined at the same time
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
  else if (strcmp(command, "sink_motion") ==0) {

    if(sink_flag == 0) error->all(FLERR,"sink_motion needs to be defined after the sink command");
    if(narg < 3) error->all(FLERR,"illegal sink_motion command");
    if(sinkmotion_flag == 0)  memory->create(sink_dr, nsink+1, "app/seg:sink_dr");
    if(sinkmotion_flag == 0)  memory->create(sink_dt, nsink+1, "app/seg:sink_dt");
    if(sinkmotion_flag == 0)  memory->create(sink_dt_new, nsink+1, "app/seg:sink_dt_new");
    if(sinkmotion_flag == 0)  memory->create(sink_dt_old, nsink+1, "app/seg:sink_dt_old");
    sinkmotion_flag = 1;

    for (i = 0; i < nsink+1; i++) {sink_dr[i] = -1.0; sink_dt[i] = 0.0;} // by default no sink motion

    int nseg = narg/3;
    for (i = 0; i < nseg; i++) {
        int sinkid = atoi(arg[i*3]); // sinkid start from 1
        if(sinkid > nsink) error->all(FLERR,"sink_motion needs to be defined after the sink command");
        sink_dr[sinkid] = atoi(arg[i*3+1]); // direction of sink motion (1 for x, 2 for y, and 3 or z), each step the displacement is a0/2
        double velocity = atof(arg[i*3+2]); // a0/s
        sink_dt[sinkid] = 0.5e12/velocity; // timestep for sink motion, in picosecond; the displacement each step is 0.5
    }
  }

  // reactions for absorption and emission
  else if (strcmp(command, "reaction") ==0) {

    if(narg != 6) error->all(FLERR,"illegal reaction command");
    if(diffusionflag) error->warning(FLERR,"MSD calculated with reactions");
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

  // calculate ris for designated elements
  else if (strcmp(command, "calris") ==0) {
    seg_flag = 1;

    if(nelement <= 0) error->all(FLERR,"calris: no elements have been defined!");
    memory->create(ris_type,nelement, "app/seg:ris_type");
    memory->create(ris_ci,nelement, "app/seg:ris_ci"); // nominal concentration as the reference
    memory->create(ris_total,nelement, "app/seg:ris_total");
    int iarg = narg/2;

    for (i = 0; i < nelement; i++) {ris_ci[i] = 0.0; ris_total[i] = 0.0;}
    for (i = 0; i < iarg; i++) {
      ris_type[i] = atoi(arg[i*2]);
      ris_ci[ris_type[i]] = atof(arg[i*2+1]);
    }
  }

  // frenkel pair production
  else if (strcmp(command, "frenkelpair") ==0) {
    frenkelpair_flag = 1;

    if(narg < 1) error->all(FLERR,"illegal frenkelpair command");

    double dose_rate = atof(arg[0]);// dose rate
    fpdistance = 0.0;
    nself_ion = 0;
    self_ion_ratio = 0.0;
    if(narg == 2) fpdistance = atof(arg[1]);// Frenkel pair separation
    if(narg == 4) {
	    self_ion_ratio = atof(arg[2]);// self-ion versus dpa
	    self_ion_type = atof(arg[3]);// self-ion type
    }
    //fpdistance *= fpdistance;

    fpfreq = 1e12/nlocal/dose_rate; // time interval to introduce an FP
    if(min_fpfreq > fpfreq) min_fpfreq = fpfreq;
  }

  // ballistic mixing
  else if (strcmp(command, "ballistic") ==0) {
    ballistic_flag = 1;

    if(narg < 2) error->all(FLERR,"illegal ballistic command");
    grow_ballistic();

    double mix_rate = atof(arg[0]);// mixing rate (in unit of dpa/s)
    bdistance = 0.0;
    bdistance = atof(arg[1]);// mixing range

    bfreq[nballistic] = 1e12/nlocal/mix_rate; // time interval for mixing
    if(min_bfreq > bfreq[nballistic]) min_bfreq = bfreq[nballistic];
    nballistic ++; // number of mixing events
  }

// meam hop steps before absorption by external sink (Peng: mfpflag)
  else if (strcmp(command, "mfp") ==0) {

    if (narg < 2 || narg % 2 == 0) error->all(FLERR,"Illegal mfp command");

    mfpflag = 1;
    memory->create(mfp,nelement+1,"app/seg:mfp");
    memory->create(nhmfp,nelement+1,"app/seg:nhmfp");
    memory->create(rhmfp,nelement+1,"app/seg:rhmfp");

    for (i = 0; i < nelement; i++) mfp[i] = -1.0;

    sigmamfp = atof(arg[0]);
    for (i=1; i<narg-1; i++) {
      if(i % 2 == 1) { j = atoi(arg[i]);
        mfp[j] = atof(arg[i+1]);
      }
    }
  }

  else error->all(FLERR,"Unrecognized command");
}

/* ----------------------------------------------------------------------
   define 2NN list based on 1NN list. In BCC lattice, each 2NN of site i is shared
   by 4 1NNs of site i as their 1NNs.
------------------------------------------------------------------------- */

void AppSeg::define_2NN()
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

void AppSeg::grow_app()
{
  type = iarray[0];   // lattice type; i1 in input
  element = iarray[1];  // element type; i2 in input
  siatype = iarray[2]; // type of sias
  dmb1 = iarray[3]; // dumbbell atom 1
  dmb2 = iarray[4]; // dumbbell atom 2

  if(diffusionflag)   aid = iarray[5]; // initially set as global ID, must use set i3 unique in command line
  if(diffusionflag || concentrationflag)   disp = darray; // msd; zero initially
}

/* ----------------------------------------------------------------------
   define virtual function site_energy(int)
------------------------------------------------------------------------- */

double AppSeg::site_energy(int)
{
  return 0.0;
}

/* ----------------------------------------------------------------------
   initialize before each run
   check validity of site values
------------------------------------------------------------------------- */

void AppSeg::init_app()
{
  int i,j;

  // constant KBT
  double KB = 0.00008617;
  KBT = temperature * KB;

  // define second nearest neighbor
  if (engstyle == 2) define_2NN();

  if (firsttime) {
    firsttime = 0;

    echeck = new int[nlocal];
    memory->create(firstevent,nlocal,"app:firstevent");

    // esites must be large enough for 3 sites and their 1st neighbors
    esites = new int[3 + 3*maxneigh];
  }

  // site statistics and sink reserved elements
  int flag = 0;
  for ( i = 0; i < nelement; i++) nsites_local[i] = 0;

  for ( i = 0; i < nlocal; i++) {
    if (type[i] < BCC || type[i] > WALL) (FLERR,"One or more sites have invalid crystal type");
    if (element[i] < VAC || element[i] > CE8) (FLERR,"One or more sites have invalid chemical type");
    nsites_local[element[i]]++;
  }

  // calculate initial global concentration
  int ntotal = 0;
  for ( i = 0; i < nelement; i++) ntotal += nsites_local[i];
  for ( i = 0; i < nelement; i++) ci[i] = 1.0*nsites_local[i]/ntotal; // need correction for V and SIA

  // initialize predefined INT by adding the two dumbbell atoms proportional to the concentrations
  dumbbell_fraction();
  set_dumbbell();
  for (i = 0; i < nelement; i++ ) {
    for (j = 0; j < nelement; j++ ) {
      emdumbbell[i][j] = 0.0;
    }
  }

  //check if reactions need to be enabled or disabled
  if(reaction_flag) {
    for( i = 0; i < nreaction; i++) {
      rcount[i] = 0;
      renable[i] = 0;
      target_global[i] = 0;
      target_local[i] = 0;
    }
    //check_reaction();
  }

  // initial recombination, elements created on sinks will be absorbed
  nFPair = 0;
  for (i = 0; i < nelement; i++) nrecombine[i] = 0;
  for (i = 0; i < nlocal; i++) recombine(i);

  // initialize the time_list for ballistic mixing
  if(ballistic_flag) {
    for(i = 0; i < nballistic; i ++) {
       time_old[i] = 0;
       time_new[i] = 0;
    }
  }

  // initialize mfp calculations(Peng: mfpflag)
  if(mfpflag) {
    for (i = 0; i < nlocal+nghost; i++) hcount[i] = 0;
    for (i = 0; i < nelement; i++) {
    //for (i = 0; i <= nelement; i++) {
        rhmfp[i] = 0.0;
        nhmfp[i] = 0;
    }

    varmfp = 2.0*sigmamfp; //sqrt(2*pi)*sigma
    sigmamfp *= sigmamfp;
    sigmamfp *= 2.0; //2*sigma^2
  }

  // initialize the time_list for sink motion
  if(sinkmotion_flag) {
    for(i = 0; i < nsink+1; i ++) {
       sink_dt_new[i] = 0.0;
       sink_dt_old[i] = 0.0;
    }
  }

  // initialize the diffusion parameters
  if(diffusionflag) {
    for(i = 0; i < 4; i ++) { // four component of diffusion
      for(j = 0; j < nlocal; j++){
         disp[i][j] = 0.0;
      }
    }

    for(i = 0; i < nelement; i ++) { // initiate the onsager coefficients
        for(j = 0; j < nelement; j++) {
	   Lij[i][j] = 0.0;
	}
    }

    for (i = 0; i < 3; i++) {
    	for (j = 0; j < nelement; j++) {
	   total_disp[i][j] = 0.0;
        }
    }
  }

 //initialize the concentration vectors
if(concentrationflag) {
  //concentration_field();
  dt_new = 0.0;
  for(i = 0; i < nelement; i ++) {
     ct[i] = 0.0;
     ct_new[i] = 0.0;

     //initiallize initial concentration
     for(j = 0; j < nlocal; j ++) {
	ct_site[i][j] = 0.0;
	disp[i+ndiff][j] = 0.0;
     }
  }
}

 // initialization for short range order
 for(i = 0; i < nelement; i++) {
    total_neighbor[i] = 0;
    for(j = 0; j < nelement; j++) {
       sro[i][j] = 0.0;
    }
 }

 // initialize solute trapping
 for(i = 0; i < nelement; i++) trap[i] = 0;

/*
  // initiate parameters for vacancy trapping
  if(time_flag) {
    double nprcs = (double)(domain->nprocs);
    fvt = 1.0/nprcs;
    itrap = 0;
    itime_current = itime_old = 0;
    treal_me = takmc_me = 0.0;
    dt_real = dt_akmc = 0.0;
  }*/

  // cacculate sink statistics
  sink_statistics();
}

/* ---------------------------------------------------------------------- */
/* ----- compute thermal sia fractions based on formation energy -------- */

void AppSeg::dumbbell_fraction()
{
 int i,j;

 double exp_total = 0.0;
 for (i = CE1; i < nelement; i ++) {
     for (j = i; j < nelement; j ++) {
	 fdumbbell[i][j] = exp(-edumbbell[i][j]/KBT);
         exp_total += fdumbbell[i][j];
     }
 }

 for (i = CE1; i < nelement; i ++) {
     for (j = i; j < nelement; j ++) {
	 fdumbbell[i][j] /= exp_total;
         fdumbbell[j][i] = fdumbbell[i][j];
     }
 }

 return;
}

/* ---------------------------------------------------------------------- */
/* ------- set initial dumbbells according to the composition ----------- */

void AppSeg::set_dumbbell()
{
 int i,j;
 double ctotal = 0.0;

 if(nsites_local[INT] == 0) return;
 for (i = CE1; i < nelement; i++) ctotal += ci[i];

 for (i = 0; i < nlocal; i++) {
     if(element[i] != INT) {
	      dmb1[i] = -1;
        dmb2[i] = -1;
	      siatype[i] = -1;
     } else {
       siatype[i] = static_cast<int>(ranseg->uniform()*3);

	     j = CE1;
       double randme1 = ranseg->uniform()*ctotal; //get 2 random # to find 2 dumbbell atoms
       double randme2 = ranseg->uniform()*ctotal;
       while(dmb1[i] <= 0) {
          if(randme1 < ci[j]) dmb1[i] = j;
          randme1 -= ci[j];
          j++;
          if(j==nelement) dmb1[i] = nelement - 1;
        }

      j = CE1;
      while(dmb2[i] <= 0) {
        if(randme2 < ci[j]) dmb2[i] = j;
        randme2 -= ci[j];
        j++;
        if(j==nelement) dmb2[i] = nelement - 1;
      }
    }
 }

 return;
}

/* ---------------------------------------------------------------------- */
/* ---------------------------initiation of looing-up list -------------- */

void AppSeg::setup_app()
{
  for (int i = 0; i < nlocal; i++) echeck[i] = 0;

  // clear event list
  nevents = 0;
  for (int i = 0; i < nlocal; i++) firstevent[i] = -1;
  for (int i = 0; i < maxevent; i++) events[i].next = i+1;
  freeevent = 0;

  // flag check and temperature
  if(temperature == 0.0) error->all(FLERR,"Temperature cannot be 0.0 for app seg");
  if(nn1flag == 0) error->all(FLERR, "First neareast neighbor bond not defined: AppSeg");
  if(nn2flag == 0 && engstyle == 2) error->all(FLERR, "Second nearest neighbor bond not defined: AppSeg");
  if(barrierflag == 0) error->warning(FLERR, "Diffusion barrier not defined: AppSeg");

  // simulation cell dimension and volume
  periodicity[0] = domain->xperiodic;
  periodicity[1] = domain->yperiodic;
  periodicity[2] = domain->zperiodic;
  lprd[0] = domain->xprd;
  lprd[1] = domain->yprd;
  lprd[2] = domain->zprd;

  boxlo[0] = domain->boxxlo;
  boxlo[1] = domain->boxylo;
  boxlo[2] = domain->boxzlo;

  volume = lprd[0]*lprd[1]*lprd[2];
}

/* ----------------------------------------------------------------------
   compute energy of site i
------------------------------------------------------------------------- */
double AppSeg::sites_energy(int i, int estyle)
{
  int j,jd,ejd,n1nn;
  double eng = 0.0;
  double cij = 0.0;

  if(estyle == 0) return eng;

  //energy from 1NN bonds
  int ei = element[i];
  n1nn = numneigh[i];  //num of 1NN
  for (j = 0; j < n1nn; j++) {
    jd = neighbor[i][j];
    ejd = element[jd];

    //adjust A-V bond based on local concentration, currently added for solute trapping
    if(trap[ei] && ejd == VAC) {
       cij = site_concentration(i,estyle);
       eng += cij*ebond1[ei][ejd];
    } else if(trap[ejd] && ei == VAC) {
       cij = site_concentration(jd,estyle);
       eng += cij*ebond1[ei][ejd];
    } else {
       eng += ebond1[ei][ejd];
    }
  }

  //energy from 2NN bonds
  if (estyle == 2) {
    int n2nn = numneigh2[i];
    for (j = 0; j < n2nn; j++) {
      jd = neighbor2[i][j];
      ejd = element[jd];

      //adjust A-V bond based on local concentration
      if(trap[ei] && ejd == VAC) {
         cij = site_concentration(i,estyle);
         eng += cij*ebond2[ei][ejd];
      } else if(trap[ejd] && ei == VAC) {
         cij = site_concentration(jd,estyle);
         eng += cij*ebond2[ei][ejd];
      } else {
         eng += ebond2[ei][ejd];
      }
    }
  }

  //bond energy shared equally by i & j
  return eng/2.0;
}

/*double AppSeg::sites_energy(int i, int estyle)
{
  int j,jd,ejd,n1nn;
  double eng = 0.0;
  double cij = 0.0;

  if(estyle == 0) return eng;

  //energy from 1NN bonds
  int ei = element[i];
  n1nn = numneigh[i];  //num of 1NN
  for (j = 0; j < n1nn; j++) {
    jd = neighbor[i][j];
    ejd = element[jd];

    //adjust A-V bond based on local concentration, currently added for solute trapping (Peng: sia interaction, need change)
    if(trap[ei] && ejd == VAC) {
       cij = site_concentration(i,estyle);
       eng += cij*ebond1[ei][ejd];
    } else if(trap[ejd] && ei == VAC) {
       cij = site_concentration(jd,estyle);
       eng += cij*ebond1[ei][ejd];
    } else if(ei == INT && ejd == INT) {
        //eng += ebond1[dmb1[i]][dmb1[jd]]+ebond1[dmb1[i]][dmb2[jd]]+ebond1[dmb2[i]][dmb1[jd]]+ebond1[dmb2[i]][dmb2[jd]];
        //total interaction: sitei dmb atom to site j dmb atom
    } else if (ejd == INT) {
      fprintf(screen, "ejd %d %d\n", dmb1[jd], dmb2[jd]);
      eng += ebond1[dmb1[jd]][ei]+ebond1[dmb2[jd]][ei];
      //total interaction: site i to dmb1 atom + site i to dmb2 atom + site i to other 1nns
    } else if (ei == INT) {
      //eng += ebond1[dmb1[i]][ejd]+ebond1[dmb2[i]][ejd]+ebond1[dmb1[i]][dmb2[i]];
      //total interaction: site i dmb1 atom to dmb2 atom + site i dmb1 to other neighbors + site i dmb2 to other 1nns
    } else {
       eng += ebond1[ei][ejd];
       // total interaction: site i to 1nns
    }
  }

  //energy from 2NN bonds
  if (estyle == 2) {
    int n2nn = numneigh2[i];
    for (j = 0; j < n2nn; j++) {
      jd = neighbor2[i][j];
      ejd = element[jd];

      //adjust A-V bond based on local concentration
      if(trap[ei] && ejd == VAC) {
         cij = site_concentration(i,estyle);
         eng += cij*ebond2[ei][ejd];
      } else if(trap[ejd] && ei == VAC) {
         cij = site_concentration(jd,estyle);
         eng += cij*ebond2[ei][ejd];
      } else {
         eng += ebond2[ei][ejd];
      }
    }
  }

  //bond energy shared equally by i & j
  return eng/2.0;
}*/

/* ----------------------------------------------------------------------
  compute local concentration to adjust solute-V bond energy
------------------------------------------------------------------------- */

double AppSeg::site_concentration(int i, int estyle)
{
  double cij = 0.0;
  int n1nn = numneigh[i]; //num of 1NN

  for(int j = 0; j < n1nn; j++) {
     int jd = neighbor[i][j];
     if(element[jd]==VAC) cij ++;
  }

  if(estyle == 1)  {
     cij = 0.5 - (cij-1)/n1nn;
     if(cij < 0.0) cij = 0.0; // attract to up to half of n1nn solutes
     return 2*cij;
  }

  int n2nn = numneigh2[i];  //num of 2NN
  for(int j = 0; j < n2nn; j++) {
     int jd = neighbor2[i][j];
     if(element[jd] == VAC) cij ++;
  }

  cij = 0.5 - (cij-1)/(n1nn+n2nn);
  if(cij < 0.0) cij = 0.0;
  return 2*cij;

}

/* ----------------------------------------------------------------------
   compute elastic energy
------------------------------------------------------------------------- */

double AppSeg::elastic_energy(int i, int itype)
{
  double eng = 0.0;
  double pressure = 0.0;

  pressure = -(stress[i][0] + stress[i][1] + stress[i][2])/3.0;
  eng = pressure * evol[itype];

  return eng;
}
/* ----------------------------------------------------------------------
  compute barriers for an exchange event between i & j except sias
------------------------------------------------------------------------- */

double AppSeg::site_SP_energy(int i, int j, int estyle)
{
  double eng = 0.0;
  double eng0i, eng0j, eng1i, eng1j; // energy before and after jump
  int ei = element[i];
  int ej = element[j];

  eng0i = sites_energy(i,estyle); //total bonds with i initially,
  eng0j = sites_energy(j,estyle); //total bonds with j initially

  // switch the element and recalculate the site energy
  element[i] = ej;
  element[j] = ei;
  eng1i = sites_energy(i,estyle); //total bonds with i after switch
  eng1j = sites_energy(j,estyle); //total bonds with j after switch

  // switch back
  element[j] = ej;
  element[i] = ei;

  //for vacancy the barrier is given by the element to be switched;
  if(ei == VAC) eng = mbarrier[ej] + eng1i + eng1j - eng0i -eng0j;

  //for other atoms the diffusion barrier is given by itself
  if(ei > INT) eng = mbarrier[ei] + eng1i + eng1j - eng0i -eng0j;

  //Contribution from segregation energy difference before and after switch; defects one step away from sink will automatically jump to sink
  if(eisink_flag && (isink[i] > 0 || isink[j] > 0)) {
    double eij = eisink[ei][isink[j]];
    double eii = eisink[ei][isink[i]];
    double eji = eisink[ej][isink[i]];
    double ejj = eisink[ej][isink[j]];

    if(eij <= -100) eij = 0.0;
    if(eii <= -100) eii = 0.0;
    if(eji <= -100) eji = 0.0;
    if(ejj <= -100) ejj = 0.0;

    eng += (eij - eii + eji - ejj)/2.0;
  }

  //add elastic contribution if applicable
  if(elastic_flag)
    eng += (elastic_energy(j,ei) - elastic_energy(i,ei) + elastic_energy(i,ej) - elastic_energy(j,ej))/2.0;

  return eng;
}

/* ----------------------------------------------------------------------
  compute barriers for an exchange event between i & j, with i being an SIA
------------------------------------------------------------------------- */

double AppSeg::sia_SP_energy(int i, int j, int estyle)
{
  double eng = 0.0;
  double eng0i, eng0j, eng1i, eng1j; //energy before and after jump
  int ei = element[i];
  int ej = element[j];
  int sia[3];
  int itype = siatype[i];
  double *vect=vdumbbell[itype];

  // check if there is an i->j diffusion path
  double dij[3];
  for (int k = 0; k < 3; k++) { // vector dij
      dij[k] = xyz[j][k] - xyz[i][k];
      if (periodicity[k] && dij[k] >= lprd[k]/2.0) dij[k] -= lprd[k];
      if (periodicity[k] && dij[k] < -lprd[k]/2.0) dij[k] += lprd[k];
  }

  //double prdt = dij[0]*vect[0] + dij[1]*vect[1] + dij[2]*vect[2];
  if (dij[itype] == 0.0) return -1.0; // dumb1 diffuses to upper or lower planes//***

  // calcualte barrier if i->j exists

  // determine which dumbbell atom diffuses
  int m = 1; // m diffuses
  int n = 2; // n stays at site i
  if(dij[itype] < 0) {m = 2; n = 1;}

  //if(itype == 0 && dij[0] < 0) m = 2;
  //if(itype == 1 && dij[1] < 0) m = 2;
  //if(itype == 2 && dij[2] < 0) m = 2;

  sia[1] = dmb1[i];
  sia[2] = dmb2[i];

  eng0i = sites_energy(i,estyle); //total bonds with i initially,
  eng0j = sites_energy(j,estyle); //total bonds with j initially
  eng0i += edumbbell[sia[1]][sia[2]]/2.0; // add sia formation energy

  // remove a dumbbell at i with the element m stays at i (m == sia1)
  element[i] = sia[n];
  dmb1[i] = -1;
  dmb2[i] = -1;
  siatype[i] = -1;

  element[j] = ei;
  eng1i = sites_energy(i,estyle); // total bonds with i after switch
  eng1j = sites_energy(j,estyle); // total bonds with j after switch
  eng1i += edumbbell[ej][sia[m]]/2.0; // dumbbell at j contains ej and sia[m]

  //fprintf(screen,"eng %f %f %f %f \n", eng0i, eng0j, eng1i, eng1j);

  // switch back
  element[j] = ej;
  element[i] = ei;

  dmb1[i] = sia[1];
  dmb2[i] = sia[2];
  siatype[i] = itype;

  //for SIA the diffusion is given by itself, correction may be added for different types of interstitials
  eng = mbarrier[ei] + eng1i + eng1j - eng0i -eng0j;
  //fprintf(screen,"eng %f %f %f %f %f \n", eng, eng0i, eng0j, eng1i, eng1j);

  //Contribution from segregation energy difference before and after switch; defects one step away from sink will automatically jump to sink
  if(eisink_flag && (isink[i] > 0 || isink[j] > 0)) {
    double eij = eisink[ei][isink[j]];
    double eii = eisink[ei][isink[i]];
    double eji = eisink[ej][isink[i]];
    double ejj = eisink[ej][isink[j]];

    if(eij <= -100) eij = 0.0;
    if(eii <= -100) eii = 0.0;
    if(eji <= -100) eji = 0.0;
    if(ejj <= -100) ejj = 0.0;

    eng += (eij - eii + eji - ejj)/2.0;
  }
  //fprintf(screen,"eng %f \n", eng);
  //add elastic contribution if applicable
  if(elastic_flag)
    eng += (elastic_energy(j,ei) - elastic_energy(i,ei) + elastic_energy(i,ej) - elastic_energy(j,ej))/2.0;
  //fprintf(screen,"eng %f \n", eng);
  //fprintf(screen,"eng %f %f %f %f %f %f\n", eng, mbarrier[ei], eng0i, eng0j, eng1i, eng1j);
  //add barrier correction by averaging starting and ending sia, in order to satisfy detailed balance
  eng += (emdumbbell[sia[1]][sia[2]] + emdumbbell[ej][sia[m]])/2;
  return eng;
}

/* ----------------------------------------------------------------------
   KMC method
   compute total propensity of owned site summed over possible events
------------------------------------------------------------------------- */

double AppSeg::site_propensity(int i)
{
  int j, k, iid, jid, sflag;

  clear_events(i);
  int ei = element[i];
  double prob_reaction = 0.0;
  double prob_hop = 0.0;
  double ebarrier = 0.0;
  double hpropensity = 0.0;
  double deltaE = 0.0;

  // Check possible reactions with barriers
  if(reaction_flag) { //reaction flag
    for(j = 0; j < nreaction; j++) {
      if(renable[j] == 0) continue;
      if(ei == rinput[j] && type[i] == rsite[j]) {
        iid = rinput[j];
        jid = routput[j];

        if(eisink_flag && eisink[element[jid]][isink[i]] != 0.0) continue; //no production at sink
       	deltaE = -sites_energy(i,engstyle); //site energy before reaction
        element[i] = jid;
        deltaE += sites_energy(i,engstyle);  //site energy after reaction
        element[i] = ei;

        ebarrier = rbarrier[j] + deltaE;
        hpropensity = rrate[j] * exp(-ebarrier/KBT);
        add_event(i,jid,2,j,hpropensity);
        prob_reaction += hpropensity;
      }
    }
  }

  if (ei > INT) return prob_reaction;

  // for hopping event propensity, the barrier is calculated by site_SP_energy();
  if (ei == VAC) { // vacancy hopping
    for (j = 0; j < numneigh[i]; j++) {
      jid = neighbor[i][j];
      if(type[jid] == WALL) continue; // no crossing the wall
      if(element[jid] <= INT) continue; // no vacancy-SIA exchange;
      ebarrier = site_SP_energy(i,jid,engstyle); // diffusion barrier

      if(ebarrier >= 0) {
        hpropensity = exp(-ebarrier/KBT);
        add_event(i,jid,1,-1,hpropensity);
        prob_hop += hpropensity;
      }
    }
  } else if (ei > VAC) {// SIA hopping
    for (j = 0; j < numneigh[i]; j++) {
      jid = neighbor[i][j];
      if(type[jid] == WALL) continue; // no crossing the wall
      if(element[jid] == INT) continue; // no SIA-SIA exchange;

      ebarrier = sia_SP_energy(i,jid,engstyle);
      if(ebarrier >= 0) {
        hpropensity = exp(-ebarrier/KBT);
      	add_event(i,jid,1,1,hpropensity);
        prob_hop += hpropensity;
      }
      //ebarrier = sia_SP_energy(i,jid,2,engstyle); // dmb2 stays at i
      //hpropensity = exp(-ebarrier/KBT);
      //add_event(i,jid,1,dmb1[i],hpropensity);
      //prob_hop += hpropensity;
    }
  }
  //fprintf(screen,"ebarrier %f %i %i\n", ebarrier, i, jid);

  return prob_hop + prob_reaction;
}

/* ----------------------------------------------------------------------
   KMC method: choose and perform an event for site i, which is selected by the solver
------------------------------------------------------------------------- */

void AppSeg::site_event(int i, class RandomPark *random)
{
  int j,k,l,m,n,ii;

  // perform events with non_zero barriers
  // pick one event from total propensity by accumulating its probability
  // compare prob to threshhold, break when reach it to select event
  double threshhold = random->uniform() * propensity[i2site[i]];
  double proball = 0.0;

  // find the event to perform
  int ievent = firstevent[i];
  while (1) {
    proball += events[ievent].propensity;
    if (proball >= threshhold) break;
    ievent = events[ievent].next;
  }

  // perform hop or reaction event
  int rstyle = events[ievent].style;
  int which = events[ievent].which; // type of reactions or neighbor id for acceleration
  j = events[ievent].jpartner;

  //fprintf(screen,"dumbbell0 %i %i\n", dmb1[i], dmb2[i]);
  //fprintf(screen,"eng0 %f %f %f\n", edumbbell[dmb1[i]][dmb2[i]], sites_energy(i,engstyle), total_energy());
  // switch element between site i and jpartner for hop diffusion
  // style 3 and 4 not currently used
  if(rstyle == 1 || rstyle == 3 || rstyle ==4) {
    k = element[i];
    if (rstyle == 4) { // switch with a 1NN of which, accelerated diffusion event
       l = element[which];
       element[i] = l;
       element[which] = element[j];
       element[j] = k;
    } else { // switch with a 1NN of i
      if(k==VAC) { //VAC mechanism
        element[i] = element[j];
        element[j] = k;
      } else { //SIA switch
	      count_dumbbell(i); // count the number of each chemical type of sias; enable when needed.
        SIA_switch(i,j);
      }

    }
    // Peng: mfpflag
    hcount[i] ++;
    if(mfpflag && mfp[element[j]] > 0.0) {
      k = hcount[i];
      hcount[i] = hcount[j];
      hcount[j] = k;

      if(hcount[j] > mfp[element[j]]) { // the defect will be absorbed
        mfp_absorption(j);
     }
   }

    // calculate MSD for each atom if activated. This is for vacancy diffusion only
    // This is valid only without defect generation and annihilation
    // MDS calcualtion for SIA diffusion is taken care of in SIA_switch()
    if(diffusionflag && element[j] == VAC) {
      // switch global atomic id
      k = aid[i];
      aid[i] = aid[j];
      aid[j] = aid[i];

      double dij[3];
      for (k = 0; k < 3; k++) { //update
        dij[k] = xyz[j][k] - xyz[i][k];
        if (periodicity[k] && dij[k] >= lprd[k]/2.0) dij[k] -= lprd[k];
        if (periodicity[k] && dij[k] <= -lprd[k]/2.0) dij[k] += lprd[k];
        disp[k][i] += dij[k];
        disp[k][j] -= dij[k];

	//count total displacement for each element for onsager calculations
	total_disp[k][element[i]] += dij[k];
	total_disp[k][element[j]] -= dij[k];
      }
      /*
      for (k = 0; k < 3; k++) { //switch
         dij[k] = disp[k][i];
         disp[k][i] = disp[k][j];
         disp[k][j] = dij[k];
      }
      disp[3][i] = disp[0][i]*disp[0][i] + disp[1][i]*disp[1][i] + disp[2][i]*disp[2][i];
      disp[3][j] = disp[0][j]*disp[0][j] + disp[1][j]*disp[1][j] + disp[2][j]*disp[2][j];
    */}

  } else {// reaction events with non-zero barriers

    k = element[i];
    element[i] = j;
    rcount[which] ++;
    nsites_local[k] --;
    nsites_local[j] ++;

    // update reaction target number
    for(ii = 0; ii < nreaction; ii++) {
      if(routput[ii] == k) target_local[ii] ++;
      if(routput[ii] == j) target_local[ii] --;
    }
  }

  //fprintf(screen,"dumbbell1 %i %i\n", dmb1[j], dmb2[j]);
  //fprintf(screen,"eng1 %f %f\n", sites_energy(i,engstyle), total_energy());
  // perform zero_barrier events: absorption and recombination
  int rid = recombine(i); // recombine site i with its neighbor rid
  if(rid >= 0) update_propensity(rid);

  //recombine j with its neighbors too after hopping;
  if(rstyle == 1 || rstyle == 3 || rstyle == 4) {
    int rid = recombine(j); // recombine site i with its neighbor rid
    if(rid >= 0) update_propensity(rid);

    //sink absorption of element[j] after hopping,for hop only since no element produced at its sinks
    if(nsink > 0) {
      absorption(i);
      absorption(j);
    }
  }

  // compute propensity changes for participating sites i & j and their neighbors
  // note that when solving using sector, j or neighbors of i may not belong to the current sector so not updated
  // recommend to solve with "sector no" option

  update_propensity(i);
  if(rstyle == 1 || rstyle == 3 || rstyle == 4) update_propensity(j);

  //LC check
  if (element[j] == element[i]){
    fprintf(screen,"ievent: %d,i:%d, element[i]:%d, jid:%d, element[j]:%d propensity:%e \n",ievent,i, element[i], j, element[j], events[ievent].propensity);
}
  // check if any reactions needs to be disabled or enabled
  if(reaction_flag == 1) {
    for(ii = 0; ii < nreaction; ii ++) {
      if(target_local[ii] <= 0 && renable[ii] == 1) {
        renable[ii] = 0;
        reset_propensity();
      }
      if(target_local[ii] > 0 && renable[ii] == 0) {
        renable[ii] = 1;
        reset_propensity();
      }
    }
  }

  return;
}

/* ----------------------------------------------------------------------
   update propensity for site i,j and their neighbors after exchange
   ignore update of sites with i2site < 0; updated during each sweep
   use echeck[] to avoid resetting propensity of same site
------------------------------------------------------------------------- */
void AppSeg::update_propensity(int i)
{
  int m,n;
  int nsites = 0;
  int isite = i2site[i];

  if (isite < 0) return;

  propensity[isite] = site_propensity(i);
  esites[nsites++] = isite;
  echeck[isite] = 1;

  for (n = 0; n < numneigh[i]; n++) { // update propensity for 1NN of site i
    m = neighbor[i][n];
    isite = i2site[m];
    if (isite >= 0 && echeck[isite] == 0) {
      propensity[isite] = site_propensity(m);
      esites[nsites++] = isite;
      echeck[isite] = 1;
    }
  }

  if(engstyle == 2) {
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

  return;
}

/* ----------------------------------------------------------------------
   check if any reaction is need to be enabled or disabled
------------------------------------------------------------------------- */
void AppSeg::check_reaction()
{
  int i,m,n,n_me,n_total,flag_local,nsites_global;
  double nprcs = (double)(domain->nprocs);

  //update the statistics of each elements locally
  for (i = 0; i < nelement; i++) nsites_local[i] = 0;
  for (i = 0; i < nlocal; i++) nsites_local[element[i]]++;

  //compute how many atoms need to be generated globally
  for (i = 0; i < nelement; i++) {
    nsites_global = 0;
    MPI_Allreduce(&nsites_local[i],&nsites_global,1,MPI_INT,MPI_SUM,world);

    for(n = 0; n < nreaction; n++){
      if(i == routput[n]) target_global[n] = rtarget[n] - nsites_global;
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
          double rand_me = ranseg->uniform();
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

  return;
}

/* ----------------------------------------------------------------------
   if any reaction is enabled or disabled, reset the propensity
   for all local lattice sites
------------------------------------------------------------------------- */
void AppSeg::reset_propensity()
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

void AppSeg::clear_events(int i)
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
   count the number of each type of SIAs for their thermodynamic fractions
------------------------------------------------------------------------- */

void AppSeg::count_dumbbell(int m)
{
  int i,j,k;

  // count by the frequencies of each type of sias (maybe more accurate by time?)
  if(siatype[m] >= 0) {
     i = ((dmb1[m]<dmb2[m]) ? dmb1[m] : dmb2[m]) - 1; //if dmb1<dmb2, i = dmb1-1; if dmb1>=dmb2, i=dmb2-1
     j = ((dmb1[m]>=dmb2[m]) ? dmb1[m] : dmb2[m]) - 1;//if dmb1>=dmb2, j = dmb1-1; if dmb1<dmb2, j=dmb2-1

     k = (i-1)*(nelement-2) + j - i*(i-1)/2 - 1;
     //nsia[k] ++; // count by frequency
     nsia[k] += dt_step; // count by residence time (*dt_step has no definition)

     /*for(i == 0; i < number_sia; i ++) {nsia[i] = 0.0;}
     for(m == 0; m < nlocal; m++) {
        if(siatype[m] >= 0) {
         i = ((dmb1[m]<dmb2[m]) ? dmb1[m] : dmb2[m]) - 1;
         j = ((dmb1[m]>=dmb2[m]) ? dmb1[m] : dmb2[m]) - 1;

            k = (i-1)*(nelement-2) + j - i*(i-1)/2 - 1;
            nsia[k] ++;
        }
     }*/

  }

  return;
}
/* ----------------------------------------------------------------------
   Perform SIA switch with an element
------------------------------------------------------------------------- */

void AppSeg::SIA_switch(int i, int j)
{
  int k;
  int ei = element[i];
  int ej = element[j];
  int itype = siatype[i];
  nsites_local[ei] --;
  nsites_local[ej] --;

  // determine the atom that diffuses
  // determine the initial and the resulting sia types and dumbbell atoms

  int jm = 1; // new dumbbellindex of m at j site
  int m = 1; // m diffuses
  int n = 2; // n stays
  int sia[3];

  sia[1] = dmb1[i];
  sia[2] = dmb2[i];

  double dij[3];
  for (k = 0; k < 3; k++) { // vector dij
      dij[k] = xyz[j][k] - xyz[i][k];
      if (periodicity[k] && dij[k] >= lprd[k]/2.0) dij[k] -= lprd[k];
      if (periodicity[k] && dij[k] < -lprd[k]/2.0) dij[k] += lprd[k];
      if (dij[k] != 0.0 && k != itype) {
	 siatype[j] = k;
	 if(dij[k] > 0) jm = 2;
      }

      // calcualte displacement for onsager coefficient calculation
      if(diffusionflag) total_disp[k][INT] += dij[k];
  }

  if(dij[itype] < 0) {m = 2; n = 1;}

  element[i] = sia[n];
  siatype[i] = -1;
  dmb1[i] = -1;
  dmb2[i] = -1;

  element[j] = ei;
  if (jm == 1) {
     dmb1[j] = sia[m];
     dmb2[j] = ej;
  } else {
     dmb2[j] = sia[m];
     dmb1[j] = ej;
  }

  // calculate displacement for onsager coefficient calculation
  // This assumes a dumbbell separation of 1/2 a0
  // The total displacement of dij[k] is splitted into three atoms
  if(diffusionflag) {
    for (k = 0; k < 3; k++) total_disp[k][sia[m]] += dij[k]/2.0; // displacement of atom that moves
    total_disp[itype][sia[n]] += dij[itype]/2.0; // atom that stays moves by half of a dumbbell separation
    total_disp[siatype[j]][ej] += dij[siatype[j]]/2.0; // the new dumbbell atom moves by half of a dumbbell separation. Displacement along the 3rd direction is zero.
  }

  //diagnose the diffusion path
  //fprintf(screen, "Debug SIA diffusion i  %d %f %f %f \n", itype,xyz[i][0],xyz[i][1],xyz[i][2]);
  //fprintf(screen, "Debug SIA diffusion j  %d %f %f %f \n", siatype[j],xyz[j][0],xyz[j][1],xyz[j][2]);

  nsites_local[element[i]] ++;
  nsites_local[element[j]] ++;

  return;
}

/* ----------------------------------------------------------------------
  add an acceleration event to list for site I
  return a propensity, Not used here (to be develeted or updated)
------------------------------------------------------------------------- */

double AppSeg::add_acceleration_event(int i, int j)
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
      if (element[nid] == VAC) { hpi[nn] = 0.0;
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
      if (element[njd] == VAC) { hpj[nn] = 0.0;
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

/* ----------------------------------------------------------------------
  add an event to list for site I
  event = exchange with site J with probability = propensity
------------------------------------------------------------------------- */

void AppSeg::add_event(int i, int j, int rstyle, int which, double propensity)
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
  //for accelerated diffusion, rstyle = 3 or 4, and switch element at site I to a neighbor of which
  //events[freeevent].kpartner = kpartner;
  events[freeevent].style = rstyle;
  events[freeevent].jpartner = j;
  events[freeevent].which = which;
  events[freeevent].propensity = propensity;

  events[freeevent].next = firstevent[i];
  firstevent[i] = freeevent;
  freeevent = next;
  nevents++;
}

/* ----------------------------------------------------------------------
   check if new Frenkal pairs need to be generated
------------------------------------------------------------------------- */
void AppSeg::check_frenkelpair(double t)
{
  int nfp = 0;
  fp_new = static_cast<int>(t/fpfreq);
  nfp = fp_new - fp_old;
  if(nfp > 100) fprintf(screen,"Too many Frenkle Pairs generated one time, %d \n", nfp);
  while (nfp > 0) {  //generating nfp frenkel pairs
    nfp --;
    frenkelpair();
    if(nfp == 0) fp_old = fp_new;  //update time
  }

  return;
}

/* ----------------------------------------------------------------------
  Create Frenkel pairs randomly. may need to scale the dose rate
  by the # of processors when work in parallel. Currently, Frenkel Pairs
  are always produced by the same processor. To be modified later.
------------------------------------------------------------------------- */
void AppSeg::frenkelpair()
{
  int i,id,vid,iid;

  // create an vacancy
  int allsites = 0;
  for (i = CE1; i < nelement; i++) allsites += nsites_local[i];
  if(allsites == 0) error->all(FLERR, "No matrix sites available for FP generation!");

  nFPair ++;
  int findv = 1;
  int velement = -1;
  while (findv) {
    id = static_cast<int> (nlocal*ranseg->uniform());
    if(id < nlocal && element[id] > INT) {
      velement = element[id];
      element[id] = VAC;
      nsites_local[VAC] ++;
      nsites_local[velement] --;

      // recalculate the propensity if defects are generated
      update_propensity(id);
      vid = id;
      findv = 0;
    }
  }

  //create an interstitial
  int findi = 1;
  int ielement = -1;
  while (findi) {
    id = static_cast<int> (nlocal*ranseg->uniform());

    if(id < nlocal && element[id] > INT) {
      if(fpdistance == 0.0) {// no requirement on fpdistance
	findi = 0;
      } else {
        double dij = distanceIJ(vid,id);
        if(dij <= fpdistance) {
	  findi = 0;
	  iid = id;
	}
      }
    }

    if(findi == 0) {
      ielement = element[iid];
      element[iid] = INT;
      dmb1[iid] = velement;
      dmb2[iid] = ielement;
      siatype[iid] = static_cast<int>(ranseg->uniform()*3);
      nsites_local[ielement] --; // site element number -1
      nsites_local[INT] ++;
      update_propensity(iid);
    }
  }

  //check if any reactions need to be activated or deactivated
  if(reaction_flag == 1) {
    for(int ii = 0; ii < nreaction; ii ++) {
      if(target_local[ii] == VAC || target_local[ii] == INT) target_local[ii] --;
      if(target_local[ii] == velement || target_local[ii] == ielement) target_local[ii] ++;

      if(target_local[ii] <= 0 && renable[ii] == 1) {
        renable[ii] = 0;
        reset_propensity();
      }
      if(target_local[ii] > 0 && renable[ii] == 0) {
        renable[ii] = 1;
        reset_propensity();
      }
    }
  }

  // no self_ions need to be added
  if(ranseg->uniform() >= self_ion_ratio) return;

  // add self-ion
  int findion = 1;
  int ionelement = -1;
  while (findion) {
    id = static_cast<int> (nlocal*ranseg->uniform());

    if(id < nlocal && element[id] > INT) {
      nself_ion ++;
      findion = 0;
      ionelement = element[id];
      element[id] = INT;
      dmb1[id] = self_ion_type;
      dmb2[id] = ionelement;
      siatype[id] = static_cast<int>(ranseg->uniform()*3);
      nsites_local[ionelement] --; // site element number -1
      nsites_local[INT] ++;
      update_propensity(id);
    }
  }

  //check if any reactions need to be activated or deactivated
  if(reaction_flag == 1) {
    for(int ii = 0; ii < nreaction; ii ++) {
      if(target_local[ii] == INT) target_local[ii] --;
      if(target_local[ii] == ionelement) target_local[ii] ++;

      if(target_local[ii] <= 0 && renable[ii] == 1) {
        renable[ii] = 0;
        reset_propensity();
      }
      if(target_local[ii] > 0 && renable[ii] == 0) {
        renable[ii] = 1;
        reset_propensity();
      }
    }
  }

  return;
}

/* ----------------------------------------------------------------------
   check if a mixing event if needed
------------------------------------------------------------------------- */
void AppSeg::check_ballistic(double t)
{
  int nmix = 0;
  for(int i = 0; i < nballistic; i ++) {
     time_new[i] = static_cast<int>(t/bfreq[i]);
     nmix = time_new[i] - time_old[i];
     while (nmix > 0) {  //perform mixing nmix times
       nmix --;
       ballistic(i);
       if(nmix == 0) time_old[i] = time_new[i];  //update time
    }
  }

  return;
}

/* ----------------------------------------------------------------------
  perform ballistic mixing. randomly mixing two atoms (no v, sia mixing)
------------------------------------------------------------------------- */

void AppSeg::ballistic(int i)
{
  int id,itype,iid,jid;

  // find atom i
  int allsites = 0;
  for(int j=CE1; j<nelement; j++) allsites += nsites_local[j];
  if(allsites == 0) error->all(FLERR, "No matrix sites available for mixing!");

  int findi = 1;
  while (findi) {
    id = static_cast<int> (nlocal*ranseg->uniform());
    if(id < nlocal && element[id] > INT) {
      iid = id;
      findi = 0;
    }
  }

  //find an element for mix
  int findj = 1;
  while (findj) {
    id = static_cast<int> (nlocal*ranseg->uniform());
    if(id < nlocal && element[id] > INT) {
        double dij = distanceIJ(iid,id);
	double mixprobability = exp(-dij/bdistance);
        if(ranseg->uniform() <= mixprobability) {
	  findj = 0;
	  jid = id;
	}
    }
  }

  //switch iid and jid and update the propensity
  itype = element[iid];
  element[iid] = element[jid];
  element[jid] = itype;
  /* currently no sia mixing
  int idmb1 = dmb1[iid];
  int idmb2 = dmb2[iid];

  dmb1[iid] = dmb1[jid];
  dmb2[iid] = dmb2[jid];

  dmb1[jid] = idmb1;
  dmb2[jid] = idmb2;
  */
  update_propensity(iid);
  update_propensity(jid);

  return;
}
/* ----------------------------------------------------------------------
   check if any sink motion needs to be performed
------------------------------------------------------------------------- */
void AppSeg::check_sinkmotion(double t)
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
------------------------------------------------------------------------- */
void AppSeg::sink_motion(int n)
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
------------------------------------------------------------------------- */
void AppSeg::absorption(int i)
{
  int j,k,m,n,ei,ejd,ntotal,rand_me,jbound,jid;
  double threshold;

  n = isink[i]; //sink id
  if(n <= 0) return; // site i is not a sink

  ntotal = 0;
  ei = element[i];
  ejd = -1;

  if(ei > INT) return; // currently only absorbs VAC and INT
  if(eisink[ei][n] > -100) return; // not an absorbing sink
  if(ranseg->uniform() > 1.0/sink_mfp[ei][n]) return; // no absorption occurs

  nabsorption[ei][n] ++; // number of element ei absorbed by type of isink[i] sink
  //fprintf(screen,"absorption %d %d \n",nabsorption[ei][n],ei);
  if(ei == VAC)  { // for vacancy choose a reserved atom to occupy
     for(m = CE1; m < nelement ; m++) {if(nreserve[m][n] > 0) ntotal += nreserve[m][n];} // count all reserved SIAs at sinks

     if(ntotal > 0) {// choose from reserved elements at sinks
       j = CE1;
       jbound = 0;
       rand_me = static_cast<int>(ranseg->uniform()*ntotal);
       while(ejd < 0) {
	 if(nreserve[j][n] > 0) jbound += nreserve[j][n];
         if(rand_me < jbound) ejd = j;
	 j++;
	 if(j == nelement) ejd = j-1;
       }
     } else {// if none reserved "borrow" one proportionally to the nominal composition from reserved
       while(ejd < 0) {
           rand_me = static_cast<int>(ranseg->uniform()*nlocal);
	   if(element[rand_me] > INT) ejd = element[rand_me];
       }
     }
     if(ejd == -1) error->all(FLERR,"no element found to recombine with vacancy at sink");

     // recombine and recount the # of atoms and reserved atomsat sink n
     nsites_local[ei] --;
     nsites_local[ejd] ++;
     element[i] = ejd;
     nreserve[ejd][n] --;

  } else if (ei == INT) { // randomly reserve one element from a dumbbell

     nsites_local[ei] --;

     // mix the two dumbbell atoms with the atoms in sink n
     // then ramdomly select one to be put into reservior
     // the probability of dumb1 or dumb2 or a sink atom to go to resevior is 1/(nsink_site+1)

     j = static_cast<int>(ranseg->uniform()*(sinksite[n-1]+2));
     if(j == sinksite[n-1]+1) { //dmb2 goes to reservior
        element[i] = dmb1[i];
        nreserve[dmb2[i]][n] ++;
     } else if(j == sinksite[n-1]) { //dmb1 goes to reservior
        element[i] = dmb2[i];
        nreserve[dmb1[i]][n] ++;
     } else { // a selected sink atom goes to reservior, dmb1 stays at i and dmb2 goes to the sink site
        element[i] = dmb1[i];
	jid = sinkid[n-1][j];
        nsites_local[element[jid]] --;
        nreserve[element[jid]][n] ++;

        element[jid] = dmb2[i];
        nsites_local[dmb2[i]] ++;
     }

     nsites_local[element[i]] ++;
     dmb1[i] = -1;
     dmb2[i] = -1;
     siatype[i] = -1;
  }

  //Peng: mfpflag reset jump steps of site i
  if(mfpflag) hcount[i] = 0;

  // update reaction target number
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

  return;
}

/*----------------------------------------------------------------------
   (Peng: mfpflag) mfp absorption of an element at site i.
-------------------------------------------------------------------------*/

void AppSeg::mfp_absorption(int i)
{
  int j,k,l,m,n,ii;

  nhmfp[element[i]] ++; //record # of defects absorbed by mean field sink
  rhmfp[element[i]] += hcount[i]; //record total # of movement for type of defect absorbed by mean field sink

  nsites_local[element[i]] --;

  //define type of element, use initial global concentration as criteria (is the criteria correct?)
  double crit = 0.0;
  //for (m = CE1; m <= nelement; m++) {
  for (m = CE1; m < nelement; m++) {
   crit += ci[m]; // here use element initial global concentration (ci[m])
   double rand_number = ranseg->uniform();
   if(rand_number < crit) {
     element[i] = m;
     break; // jump out of the for loop
   }
  }

  nsites_local[m] ++;
  hcount[i] = 0;
}

/*-----------------------------------------------------------------
recombination within 2nd NN distance, returns the id of the recombined atom with i
-------------------------------------------------------------------*/
int AppSeg::recombine(int i)
{
  int n1nn,n2nn,jd,ei,ej,m,n,sia1,sia2;

  n1nn = numneigh[i];
  n2nn = 0;
  ei = element[i];
  if(engstyle == 2) n2nn = numneigh2[i];
  for (n = 0; n < n1nn + n2nn; n++) {
      if(n < numneigh[i])  jd = neighbor[i][n];
      if(n >= numneigh[i])  jd = neighbor2[i][n-n1nn]; // 2NN

      ej = element[jd];
      if(ei != VAC && ej != VAC) continue; // None vacancy
      if(ei != INT && ej != INT) continue; // None SIA
      if(ei == ej) continue; // same type

      nrecombine[ei] ++;
      nrecombine[ej] ++;
      nsites_local[ei] --;
      nsites_local[ej] --;

      // extract the two dumbbell atoms
      if(ei == INT) {
        sia1 = dmb1[i];
        sia2 = dmb2[i];
	siatype[i] = -1;
      } else {
        sia1 = dmb1[jd];
        sia2 = dmb2[jd];
	siatype[jd] = -1;
      }

      // randomly assign the two dumbbell atoms to sites i and jd
      double randme = ranseg->uniform();
      if(randme < 0.5) {
	element[i] = sia1;
	element[jd] = sia2;
      } else {
	element[i] = sia2;
	element[jd] = sia1;
      }

      nsites_local[element[i]] ++;
      nsites_local[element[jd]] ++;

     // if(mfpflag) {hcount[i] = 0; hcount[m] = 0;}
     if(mfpflag) {hcount[i] = 0; hcount[jd] = 0;} //Peng: mfpflag

     // update reaction target number
     if(reaction_flag == 1) {
       for(int k = 0; k < nreaction; k++) {
       if(routput[k] == element[i] || routput[k] == element[jd])  target_local[k] --;
       if(routput[k] == ei || routput[k] == ej)  target_local[k] ++;
       }
     }

     return jd;
  }

  // no recombination occurs
  return -1;
}

/* ----------------------------------------------------------------------
  calculating total energy
------------------------------------------------------------------------- */

double AppSeg::total_energy( )
{
  int i,j,etype,stype;
  double penergy = 0.0;
  for(i = 0; i < nlocal; i++) penergy += sites_energy(i,engstyle);
  //diagnose if any dumbbells are assigned to non-interstitial atoms
  //for(i = 0; i < nlocal; i++) {if(element[i] != INT && (dmb1[i] >= 0 || dmb2[i] >= 0))
  //   {fprintf(screen,"wrong atoms in non-SIAs %d %d %d %d\n", i,element[i],dmb1[i],dmb2[i]);}}

  for(i = 0; i < nlocal; i++) {
     if(siatype[i] >= 0) penergy += edumbbell[dmb1[i]][dmb2[i]];
  }

  if(elastic_flag) { // elastic energy
    for(j = 0; j < nlocal; j++) {
      etype = element[j];
      penergy += elastic_energy(j,etype);
    }
  }

  if(eisink_flag) { // segregation energy
    for(j = 0; j < nlocal; j++) {
      etype = element[j];
      stype = isink[j];
      if(stype > 0 && eisink[etype][stype] > -100) penergy += eisink[etype][stype];;
    }
  }

//LC check
  // for (int i = 0; i < nevents; i++){
  //   fprintf(screen, "index: %d ,rstyle:%d, jid: %d , element[j]: %d \n",i ,events[i].style,  events[i].jpartner, element[events[i].jpartner]);
  //      }

  return penergy;
}

/* ----------------------------------------------------------------------
  calculate the Onsager coefficient based on atomic displacement
  Lij=<DR_i>*<Dr_j>/6VKbTt;
------------------------------------------------------------------------- */
void AppSeg::onsager(double t)
{
  int i,j;
  //double total_disp[3][10];

  if(t <= 0) return;
/*
  for (i = 0; i < 3; i++) {
  for (j = 0; j < nelement; j++) {
	  total_disp[i][j] = 0.0;
  }
  }

  // calculate the total displacement of each element
  for (i = 0; i < nlocal; i++) {
     total_disp[0][element[i]] += disp[0][i];
     total_disp[1][element[i]] += disp[1][i];
     total_disp[2][element[i]] += disp[2][i];
  }
*/
  for (i = 0; i < nelement; i++) {
  for (j = i; j < nelement; j++) {
      Lij[i][j] = total_disp[0][i]*total_disp[0][j] + total_disp[1][i]*total_disp[1][j]+total_disp[2][i]*total_disp[2][j];
      Lij[i][j] /=(6*t*KBT*volume*1e-12);
      Lij[j][i] = Lij[i][j];
  }
  }

  return;
}

/* ----------------------------------------------------------------------
  calculate the short range order matrix
------------------------------------------------------------------------- */
void AppSeg::short_range_order()
{ int i,j,jd,itype,jtype;

  // initialization
  for(i=0; i<nelement; i++) {
     total_neighbor[i] = 0;
     for(j=0; j<nelement; j++) {
        sro[i][j] = 0.0;
     }
  }

  for(i=0; i<nlocal; i++) {
     itype = element[i];
     for(j=0; j<numneigh[i]; j++) {
	jtype=element[neighbor[i][j]];
        total_neighbor[itype] ++;
	sro[itype][jtype] += 1;
     }
  }

  for(i=0; i<nelement; i++) {
     for(j=0; j<nelement; j++) {
        double jconcentration = 1.0*nsites_local[j]/nlocal;
	if(jconcentration == 0.0) { // there is no j element
	  sro[i][j] = 1.0;
	} else {
          sro[i][j] /= total_neighbor[i];
  	  sro[i][j] /= jconcentration;
          sro[i][j] = 1.0 - sro[i][j];
	}
     }
  }

  return;
}

/* ----------------------------------------------------------------------
  Calculate time dependent ris of given element
------------------------------------------------------------------------- */
void AppSeg::ris_time()
{
  int i,j,ncell;
  double **icell;

  ncell = static_cast<int>(lprd[2]); // bin size = 1
  icell = new double *[ncell];
  for(i = 0; i < ncell; i++) {
     icell[i] = new double[nelement];
  }

  for(i = 0; i < nelement; i++) ris_total[i] = 0.0;
  for(i = 0; i < ncell; i++){
  for(j = 0; j < nelement; j++) {
     icell[i][j] = 0;
  }}

  for(i = 0; i < nlocal; i++) {
     int iz = static_cast<int>(xyz[i][2] - boxlo[2]);
     if(concentrationflag) {
       icell[iz][element[i]] += disp[ndiff+element[i]][i];
     } else {
       icell[iz][element[i]] += 1.0;
     }
  }

  int nlayer = nlocal/ncell;
  for(j = 0; j < nelement; j++) {
     if(ris_ci[j] == 0.0) continue; // do not calcualte for those not monitored
     for(i = 0; i < ncell; i++) {
        double dcij = 1.0*icell[i][j]/nlayer - ris_ci[j]; // deviation from nominal concentration
        ris_total[j] += fabs(dcij)/2.0; // count both enrichment and depletion and the divide the sum by 2
     }
  }

  for(i = 0; i < ncell; i++) {
     delete [] icell[i];
  }
  delete [] icell;

  return;
}

/* ----------------------------------------------------------------------
  Integrate c*t at each site for fractional occupancy over time
------------------------------------------------------------------------- */
void AppSeg::concentration_field(double dt)
{
  dt_new += dt; // update time interval
  for(int i = 0; i < nlocal; i++) {
     ct_new[element[i]] += dt; // total concentration
     ct_site[element[i]][i] += dt; // site concentration
  }

  return;
}

/* ----------------------------------------------------------------------
  update time averaged total concentration concentrations
------------------------------------------------------------------------- */
void AppSeg::time_averaged_concentration()
{
  if(dt_new <= 0) return;

  for(int i = 0; i < nelement; i++) { // ct = c*t / dt
     ct[i] = ct_new[i]/dt_new/nlocal; //time and spatial average
     ct_new[i] = 0.0; //start recounting

     for(int j = 0; j < nlocal; j++) {
	disp[ndiff+i][j] = ct_site[i][j]/dt_new; //time average
	ct_site[i][j] = 0.0; //start recounting

     }
  }

  dt_new = 0.0;
}

/* ----------------------------------------------------------------------
  check if the vacancy is trapped by solute by counting its
  first and second NNs, return 0 if any of those is non-Fe
------------------------------------------------------------------------- */
int AppSeg::vacancy_trap(int i)
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
  Peng: check if the neighbor of certain B atom are all A atoms or not and count
The number of B monomers and dimers
------------------------------------------------------------------------- */
int AppSeg::solubility(int n)
{
  int j,jd,k,z;

  //check neighbor for B atoms
  di_mono = 0;
  for (j = 0; j < nlocal; j++) {
    z = 0;
    if (element[j] == n) {
      for (k = 0; k < numneigh[j];k++){
        jd = neighbor[j][k];
        if (element[jd] != n) {
        }
        else {
          z = z+1;
        }
      }

      if (z <= 1){
        di_mono = di_mono+1;
        }
    }
  }
  return di_mono;
}

/* ----------------------------------------------------------------------
  record KMC time and real time
------------------------------------------------------------------------- */
void AppSeg::time_tracer(double dt)
{
  dt_akmc += dt;
  dt_real += dt*itrap; // not contribute to real time if trapped by solute
}

/* ----------------------------------------------------------------------
  check if the vacancy is trapped by solute, contribute 0 to physical
  time if so, return zero if itrap == 1; 1 otherwise
------------------------------------------------------------------------- */

double AppSeg::real_time(double t)
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

int AppSeg::ibonde(int a, int b, int c)
{
  return ((a-1)*c + b - a*(a-1)/2);
}

/* ----------------------------------------------------------------------
  grow memory for dislocation
------------------------------------------------------------------------- */

void AppSeg::grow_dislocations()
{
  int n = ndislocation + 1;
  memory->grow(dislocation_type,n,"app/ris:dislocation_type");
  memory->grow(burgers,n,3,"app/ris:burgers");
  memory->grow(xdislocation,n,3,"app/ris:xdislocation");
  memory->grow(line_vector,n,"app/ris:line_vector");
  memory->grow(dislocation_radius,n,"app/ris:dislocation_radius");
  memory->grow(nsegment,n,"app/ris:nsegment");
}

/* ----------------------------------------------------------------------
  grow memory for sink
------------------------------------------------------------------------- */

void AppSeg::grow_sinks()
{
  int n = nsink + 1;
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

void AppSeg::grow_reactions()
{
  int n = nreaction + 1;
  memory->grow(rsite,n,"app/ris:rsite");
  memory->grow(rinput,n,"app/ris:rinput");
  memory->grow(routput,n,"app/ris:routput");
  memory->grow(rcount,n,"app/ris:rcount");
  memory->grow(renable,n,"app/ris:renable");
  memory->grow(rtarget,n,"app/ris:rtarget");
  memory->grow(rbarrier,n,"app/ris:rbarrier");
  memory->grow(rrate,n,"app/ris:rrate");
  memory->grow(target_local,n,"app/ris:target_local");
  memory->grow(target_global,n,"app/ris:target_global");
}

/* ----------------------------------------------------------------------
  grow memory for ballistic mixing
------------------------------------------------------------------------- */

void AppSeg::grow_ballistic()
{
  int n = nballistic + 1;
  memory->grow(bfreq,n,"app/ris:bfreq");
  memory->grow(time_old,n,"app/ris:time_old");
  memory->grow(time_new,n,"app/ris:time_new");
}

/* ----------------------------------------------------------------------
  create sinks
------------------------------------------------------------------------- */

void AppSeg::sink_creation(int n)
{
  int i,j,periodicity[3];
  int shape = sink_shape[n];
  int normal = sink_normal[n];
  int segment = sink_segment[n];
  int nlattice = nlocal + nghost;
  double dx,dij[3],rik,rjk,lprd[3];
  double radius = sink_radius[n];
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
      if(dx < range) isink[i] = n;
    }
  }

  return;
}

/* ----------------------------------------------------------------------
  sink_statistics, record the # sites and their ids for each sink
------------------------------------------------------------------------- */

void AppSeg::sink_statistics()
{ int i,j,maxsite,id;

  // Record the sink site ids

  maxsite = 0;
  memory->create(sinksite,nsink,"app/seg:sinksite"); // # of sink sites for each sink
  for(i = 0; i < nsink; i++) {
     sinksite[i] = 0;
     for(j = 0; j < nlocal; j++) {
        if(isink[j] == i+1) sinksite[i] ++; // sink id starts from 1
     }
     if(sinksite[i] > maxsite) maxsite = sinksite[i];
  }

  memory->create(sinkid,nsink,maxsite,"app/seg:sinkid"); // ids of each sink site for each sink
  for(i = 0; i < nsink; i++) {
     id = 0;
     for(j = 0; j < nlocal; j++) {
        if(isink[j] == i+1) {
          sinkid[i][id] = j;
	  id++;
	}
     }
  }

  return;
}
/* ----------------------------------------------------------------------
  calculate dislocation stess field
------------------------------------------------------------------------- */

void AppSeg::stress_field(int n)
{
  if(dislocation_type[n] == 1) stress_dislocation(n); //straight dislocation
  else stress_loop(n); //loop
}

/* ----------------------------------------------------------------------
  calculate stess field for straight dislocations
------------------------------------------------------------------------- */

void AppSeg::stress_dislocation(int n)
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

/* ----------------------------------------------------------------------
  calculate stess field for dislocation loops
------------------------------------------------------------------------- */

void AppSeg::stress_loop(int n)
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
/* ----------------------------------------------------------------------
  Calculate stress field due to segment AB
------------------------------------------------------------------------- */

void AppSeg::seg_stress( double A[3], double B[3], double P[3], double bv[3], double sstres[3][3])
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

/* ----------------------------------------------------------------------
  Calculate sigma_A
------------------------------------------------------------------------- */

void AppSeg::sigma_A(double t[3], double m[3], double bv[3], double sigma[3][3])
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

/* ----------------------------------------------------------------------
  Calculate sigma_P
------------------------------------------------------------------------- */

void AppSeg::sigma_P(double t[3], double ni[3], double bv[3], double sigmap[3][3])
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

/* ----------------------------------------------------------------------
  stress calculation based on stroh formulism
------------------------------------------------------------------------- */

void AppSeg::stroh( double tvect[3], double qmatx[3][3], double bmatx[3][3], double smatx[3][3])
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

/* ----------------------------------------------------------------------
  angular factors in stress calculation based on stroh formulism
------------------------------------------------------------------------- */

void AppSeg::stroh_p( double t[3], double n0[3], double qpmat[3][3], double bpmat[3][3], double spmat[3][3])
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

/* ----------------------------------------------------------------------
  establish elastic tensor matrix from c11 c12 & c44
------------------------------------------------------------------------- */

void AppSeg::elastic_tensor()
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

/* ----------------------------------------------------------------------
  3x3 matrix inversion
------------------------------------------------------------------------- */

void AppSeg::matrix_inversion(double mm[3][3], double nn[3][3])
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

void AppSeg::cross_product( double m[3], double n[3], double l[3])
{
  l[0] = m[1]*n[2] - m[2]*n[1];
  l[1] = m[2]*n[0] - m[0]*n[2];
  l[2] = m[0]*n[1] - m[1]*n[0];

}

/* ----------------------------------------------------------------------
  cross product of int vectors
------------------------------------------------------------------------- */

void AppSeg::cross_product( int m[3], int n[3], int l[3])
{
  l[0] = m[1]*n[2] - m[2]*n[1];
  l[1] = m[2]*n[0] - m[0]*n[2];
  l[2] = m[0]*n[1] - m[1]*n[0];

}

/* ----------------------------------------------------------------------
  normalize double vector
------------------------------------------------------------------------- */

void AppSeg::vector_normalize( double m[3])
{
  double norm = sqrt(m[0]*m[0] + m[1]*m[1] + m[2]*m[2]);
  if(norm == 0) error->all(FLERR,"can not normalize zero vector");
  for(int i = 0; i < 3; i ++) m[i] /= norm;

}

/* ----------------------------------------------------------------------
  right-hand coordination system, m x n x t
------------------------------------------------------------------------- */

void AppSeg::right_hand_coord( double m[3], double n[3], double t[3])
{
  double taux[3];
  // random vector taux to find m
  taux[0] = ranseg->uniform();
  taux[1] = ranseg->uniform();
  taux[2] = ranseg->uniform();

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

void AppSeg::trapozidal(double omegas[100], double XP[9][100], int j, int k, double value)
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
   calculate the distance between i and j in vector form.
------------------------------------------------------------------------- */
double AppSeg::distanceIJ(int i, int j)
{
   //update and switch displacement
   int periodicity[3];
   double lprd[3],dij[4];

   periodicity[0] = domain->xperiodic;
   periodicity[1] = domain->yperiodic;
   periodicity[2] = domain->zperiodic;
   lprd[0] = domain->xprd;
   lprd[1] = domain->yprd;
   lprd[2] = domain->zprd;

   dij[3] = 0.0;
   for (int k = 0; k < 3; k++) { //update
       dij[k] = xyz[j][k] - xyz[i][k];
       if (periodicity[k] && dij[k] >= lprd[k]/2.0) dij[k] -= lprd[k];
       if (periodicity[k] && dij[k] <= -lprd[k]/2.0) dij[k] += lprd[k];
       dij[3] += dij[k]*dij[k];
   }

   return sqrt(dij[3]); //square root of distance
}
