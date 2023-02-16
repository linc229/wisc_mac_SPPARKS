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
AppStyle(coros,Appcoros)
#else

#ifndef SPK_APP_coros_H
#define SPK_APP_coros_H

#include "app_lattice.h"

namespace SPPARKS_NS {

class Appcoros : public AppLattice {
  friend class Diagcoros;

 public:
  Appcoros(class SPPARKS *, int, char **);
  ~Appcoros();
  void input_app(char *, int, char **);
  void grow_app();
  void init_app();
  void setup_app();

  double site_energy(int);
  void site_event_rejection(int, class RandomPark *) {}
  double site_propensity(int);
  void site_event(int, class RandomPark *);

 private:
  int engstyle,nn1flag,nn2flag,barrierflag,diffusionflag; // 1NN or 2NN bonds
  int attemptfrequencyflag; // flag for attempt frequency by LC
  //int ndiffusion;
  int ndiff;
  int *type,*element,*aid; // variables on each lattice site
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
  //double *surfbarrier; //migration barrier for surface diffusion by LC comment 0704
  double *hcount;  // LC change to big int
  double *attemptfrequency; // attempt freqency for all diffusion by LC

  int nreact; // number of reaction events by LC
  //int nbulkfe; // number of bulk diffusion event for id2 = 1 by LC
  //int nbulkcu; // number of bulk diffusion event for id2 = 3 by LC

  int nbulk;
  int ninterface;
  int nsalt;

  int nrecombine; //number of recombination event

  class RandomPark *rancoros; //random number generator
//parameter for dislocaitons
  int moduli_flag,dislocation_flag,elastic_flag,ndislocation,ninteg;
  int *dislocation_type,*line_vector,*nsegment;
  double **burgers,**xdislocation,**stress;
  double *dislocation_radius;
  double c11,c12,c44,dcore;
  double evol[10],cijkl[3][3][3][3];

//oarameter for sinks // LC comments
  // int nsink,sink_flag,eisink_flag;
  // int *sink_type,*sink_shape,*sink_segment,*sink_normal,*nabsorption;
  // int **isink;
  // double *sink_strength,*sink_radius,*sink_mfp;
  // double **xsink;

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

//paramter for vacancy trapping
  int itrap,itime_current,itime_old;
  double dt_interval,fvt,dt_real,dt_akmc,treal_me,takmc_me;

//parameter for ballistic mixing
  //int saltdiffusion_flag; // LC
  int nballistic;
  int *time_old,*time_new;
  double *rdamp,*pn_local,*pn_global,*bfreq;
  double **xmix,**pmix;

//parameter for salt potential by LC
  int nsaltdiffusion; // LC
  int *potential; //atomic displacement
  bigint *salt_time_old,*salt_time_new;
  double *salt_bfreq;
  bigint num_saltdiffusion ; // number of salt diffusion motion
  int temp_salt_time; // LC

//parameter for acceleration
  int ntrap;
  int *trap_type;

  //parameter for update_region
  int n_update_list;
  int update_list[200]; //set array size 200

  //parameter for barrier_extract and data_extract
  int extract_flag;
  int evap_extract_flag;
  int np_extract_flag;
  //int coros_stop_flag;
  int threshold_Cr;
  int large_propensity_extract_flag;

  //parameter for surface_effect
  double surface_effect_b;

  //parameter for initial particle count;
  int total_Ni;
  int total_vac;
  int total_Cr;
  int surface_Ni;
  int surface_Cr;

  //parameter for time averaged concentration
    double **ct_site, **ct_site_new, *ct,*ct_new,dt_new, dt_site_c_new;
    double *ct_site_temp_t;  //LC for stoing time point
    int *ct_site_temp_i2;    //LC for storing i2 at previous time point
    int cur_i;       // global var for storing position of i and j
    int cur_j;
    //double site_time_interval; //LC
    double time_sector; // LC
    double ct_time; // LC
    double active_vac_new, active_vac; // LC
    double *site_time;                  // LC
    double *site_time_interval;        // LC
    //int ct_reset_flag;
    //int ct_site_flag;
    double **i3_site, **i3_site_new, dt_i3_site_new;
  // parameter for dump_event func by LC
  //int dump_event_flag;
  int dump_index;
  int dump_event_all_flag;
  int dump_event_all_file;
  int dump_event_all_index;

  // parameters for corrosion_info_output by LC;
  //int coros_flag;
  int coros_index;
  int *coros_id;
  double *coros_time;


  //parameter for monomer count by LC
  int *monomers;

  //parameters for msd calculation by LC
  double *sd;

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


  int ibonde(int, int, int);  //list to matrix bond energy
  void define_2NN();
  void clear_events(int);
  void add_event(int, int, int, int, double);
  void update_propensity(int);
  double add_acceleration_event(int, int); // add acceleration event and return a probability
  double total_energy();
  double sites_energy(int, int);
  double site_SP_energy(int, int, int, double);

  void grow_reactions(); //reactions
  void check_reaction();
  void reset_propensity();

  void grow_ballistic();// ballistic mixing
  void check_ballistic(double);
  void ballistic(int);
  void ballistic_probability(int);
  void absorption(int);

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

  void grow_sinks(); //sink
  void sink_creation(int);
  void check_sinkmotion(double);
  void sink_motion(int);
  void sink_statistics();

  int vacancy_trap(int);
  void time_tracer(double); //track time
  //void concentration_field(); //calculation concentration field
  double real_time(double); //compute fvt
  //void update_region(int i,int j, int r); // // update type after events
  //int update_neighbor_check(int l); //update and return number of old list
  //int update_surface_diff(int i); // update surface diff
  void count_type();// to count each type
  void potential_diff(int nsalt); //perform salt potential_diffusion
  void check_saltdiffusion(double); // check salt diffusion and time
  void grow_saltdiffusion();// grow memory for salt diffusion
  void salt_remove(int i); //remove salt potential after reaction by LC
  int count_salt(); // count salt if i3 =1
  void barrier_print(int r,double i, double j, double k, double l); //by LC
  void np_check(int i, int jid); //by LC
  double total_metal_energy(); // by LC
  int KMC_stop(); //by LC
  int vac_monomer_count(); //by LC
  void concentration_field(double); //calculation concentration field
  void concentration_field_global(double); //LC for global
  void time_averaged_concentration(); // calculate time-averaged concentration
  void site_concentration_calc(int, int);     // by LC
  void monomer_count();
  //void ct_reset(); // reset time-average concentration after diag call
  double **ct_site_extract(); // return ct_site **array called by app_lattice <-- dump_text
  double **i3_site_extract(); // return i3_site **array called by app_lattice <-- dump_text
  double *ct_extract();       // LC
  double active_vac_extract(); // LC
  int metal_pure_vac_approxi(); //LC, return pure_vac in metal region by approximation
  void dump_event(double); // LC, dump event list per dump time as text file
  void dump_event_all(int, int, int, int, double); // LC dump all event list to one file
  void MSD_calculation(); //by LC, compute MSD when atom change
  double *temp_t_extract(); // by LC
  int *temp_i2_extract();   // by LC

  void corrosion_info_output(int i);   // by LC
  void corrosion_info_dump();          // by LC
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

E: Temperature cannot be 0.0 for app coros

UNDOCUMENTED

*/
