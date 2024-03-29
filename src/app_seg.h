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

#ifdef APP_CLASS
AppStyle(seg,AppSeg)
#else

#ifndef SPK_APP_SEG_H
#define SPK_APP_SEG_H

#include "app_lattice.h"

namespace SPPARKS_NS {

class AppSeg : public AppLattice {
  friend class DiagSeg;

 public:
  AppSeg(class SPPARKS *, int, char **);
  ~AppSeg();
  void input_app(char *, int, char **);
  void grow_app();
  void init_app();
  void setup_app();

  double site_energy(int);
  void site_event_rejection(int, class RandomPark *) {}
  double site_propensity(int);
  void site_event(int, class RandomPark *);

 private:
  int engstyle,nn1flag,nn2flag,barrierflag,seg_flag; // 1NN or 2NN bonds
  int mfpflag; //Peng: mfpflag
  int ndiff,number_sia,ndumbbell;
  int *type,*element,*dmb1,*dmb2,*aid,*siatype; // variables on each lattice site
  int firsttime;

  int *esites;
  int *echeck;

  int *numneigh2; // number of 2NN
  int **neighbor2; // index of 2NN

  int nelement; //total reactions
  int *nsites_local; //statics of local element
  double KBT; //rate and propensity for hop diffusion
  double **ebond1,**ebond2; //bond energy
  double **disp; //atomic displacement
  double *mbarrier; //migration barriers
  double **vdumbbell; //dumbbell vector
  double **edumbbell; //dumbbell formation energy
  double **emdumbbell; //dumbbell migration barrier correction
  double **fdumbbell; //relative fraction of dumbbells
  double *nsia; //number of each type of dumbbells
  int *hcount;
  int periodicity[3]; //periodicity;
  double total_disp[3][10],boxlo[3],lprd[3],volume; //simulation cell size

  /*-----------------------------Peng: mfpflag-------------------------------*/
  double sigmamfp,varmfp; //rate and propensity for hop diffusion
  double *mfp,*rhmfp; //migration barriers
  int *nhmfp;
  /*=========================================================================*/

  int nrecombine[9]; //number of recombination event

  class RandomPark *ranseg; //random number generator
//parameter for dislocaitons
  int moduli_flag,dislocation_flag,elastic_flag,ndislocation,ninteg;
  int *dislocation_type,*line_vector,*nsegment;
  double **burgers,**xdislocation,**stress;
  double *dislocation_radius;
  double c11,c12,c44,dcore;
  double evol[11],cijkl[3][3][3][3];

//parameter for sinks
  int nsink,sink_flag,eisink_flag;
  int *isink,*sink_shape,*sink_segment,*sink_normal,*sinksite;
  int **nabsorption,**nreserve,**sinkid;
  double *ci,*sink_range,*sink_radius,*sink_dr,*sink_dt,*sink_dt_new,*sink_dt_old;
  double **xsink,**eisink,**sink_mfp;

//parameter for reaction
  int nreaction;
  int *rsite,*rinput,*routput,*rcount,*renable,*rtarget;
  int *target_local,*target_global;
  double *rbarrier,*rrate; // rrate is scaled by the attempt rate of atom hopping

//parameters for solute trapping to adjust bond energies
  int trap[10];

//paramter for vacancy trapping for time rescaling
  int itrap,itime_current,itime_old;
  double dt_interval,fvt,dt_real,dt_akmc,treal_me,takmc_me;

//parameter for ballistic mixing
  int nballistic, nFPair;
  int *time_old,*time_new;
  double *bfreq;
  double bdistance;

//parameter for frenkel pair generation
  int nfp,fp_old,fp_new,nself_ion,self_ion_type;
  double fpfreq,fpdistance,self_ion_ratio;

//parameter for time averaged concentration
  double **ct_site,*ct,*ct_new,dt_new;

//parameter for ris calculation
  int *ris_type;
  double *ris_total;
  double *ris_ci;

//parameter for acceleration
  int ntrap;
  int *trap_type;

//parameter for Onsager coefficient calculation
  double **Lij;

//parameter for Onsager coefficient calculation
  int *total_neighbor;
  double **sro;

  struct Event {           // one event for an owned site
    int style;             // reaction style = HOP,RECOMBINE
    int which;             // which reaction of this type
    int jpartner;          // which J neighbors of I are part of event
    int next;              // index of next event for this site
    double propensity;     // propensity of this event
  };

  Event *events;           // list of events for all owned sites
  int nevents;             // # of events for all owned sites
  int maxevent;            // max # of events list can hold
  int *firstevent;         // index of 1st event for each owned site
  int freeevent;           // index of 1st unused event in list
  int di_mono;                // Peng: number of B monomers

  int ibonde(int, int, int);  //list to matrix bond energy
  void define_2NN();
  void clear_events(int);
  void add_event(int, int, int, int, double);
  void update_propensity(int);
  double add_acceleration_event(int, int); // add acceleration event and return a probability
  double total_energy();
  double sites_energy(int, int);
  double site_concentration(int, int);
  double site_SP_energy(int, int, int);
  double sia_SP_energy(int, int, int);

  void grow_reactions(); //reactions
  void check_reaction();
  void reset_propensity();

  void grow_ballistic();// ballistic mixing
  void check_ballistic(double);
  void check_frenkelpair(double);
  void ballistic(int);
  void frenkelpair();
  void absorption(int);
  void ballistic_probability(int);
  void mfp_absorption(int);//Peng: mfpflag

  void grow_dislocations(); //dislocation
  void stress_field(int);
  void stress_loop(int);
  void elastic_tensor();
  void stress_dislocation(int);
  void vector_normalize(double []);
  void cross_product(int [], int [], int []);
  void cross_product(double [], double [], double []);
  void matrix_inversion(double [][3], double[][3]);
  void right_hand_coord(double [], double [], double []);
  void stroh(double [], double [][3], double [][3], double [][3]);
  void stroh_p(double [], double [], double [][3], double[][3], double [][3]);
  void trapozidal(double [], double[][100], int, int, double);
  void sigma_A(double [], double [], double [], double [][3]);
  void sigma_P(double [], double [], double [], double [][3]);
  void seg_stress(double [], double [], double [], double [], double [][3]);
  double elastic_energy(int,int);
  double distanceIJ(int,int); //distance between site I&J
  void dumbbell_fraction(); //thermodynamic participation of dumbbells

  void grow_sinks(); //sink
  void sink_creation(int);
  void check_sinkmotion(double);
  void sink_motion(int);
  void sink_statistics();

  void set_dumbbell(); // set the dumbbell atoms
  void count_dumbbell(int); // count the chemicaltype of each dumbbells

  int recombine(int); // recombination
  int vacancy_trap(int);
  void time_tracer(double); //time tracer
  void SIA_switch(int,int); //SIA_switch
  void concentration_field(double); //calculation concentration field
  void time_averaged_concentration(); // calculate time-averaged concentration
  double real_time(double); //compute fvt

  void ris_time(); // calculate ris
  void onsager(double); // calculate the Onsager coefficients
  void short_range_order(); // calculate the short range order
  int solubility(int); // Peng:function of B_monomer & dimer
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running SPPARKS to see the offending
line.

E: Unrecognized command

The command is assumed to be application specific, but is not
known to SPPARKS.  Check the input script.

E: One or more sites have invalid values

The application only allows sites to be initialized with specific
values.

E: Temperature cannot be 0.0 for app ris

UNDOCUMENTED

*/
