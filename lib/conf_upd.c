#ifndef GAUGE_CONF_UPD_C
#define GAUGE_CONF_UPD_C

#include"../include/macro.h"

#include<math.h>
#include<stdlib.h>

#include"../include/conf.h"
#include"../include/gparam.h"
#include"../include/random.h"

// staples for the plaquette component of the action
// sum (plaq) = lambda_{x,mu}*staple + independent of lambda_{x,mu}
int plaqstaples_for_link(Conf *GC,
                         Geometry const * const geo,
                         long r,
                         int i)
  {
  int j, k;
  int ris, tmp;
  long r1;

  ris=0;

  for(j=i+1; j<i+STDIM; j++)
     {
     k=j%STDIM;

//              ^ i
//         (6)  |  (3)
//      +---<---+--->---+
//      |       |       |
//   (5)V       ^       V (2)
//      |       |       |
//      +--->---+---<---+---> k
//     r1  (4)  r  (1)

     tmp=GC->lambda[r][k]; //(1)
     tmp*=GC->lambda[nnp(geo, r, k)][i]; //(2)
     tmp*=GC->lambda[nnp(geo,r,i)][k]; //(3)
     ris+=tmp;

     r1=nnm(geo, r, k);
     tmp= GC->lambda[r1][k]; //(4)
     tmp*=GC->lambda[r1][i]; //(5)
     tmp*=GC->lambda[nnp(geo, r1, i)][k]; //(6)
     ris+=tmp;
     }

    return ris;
    }


// perform an update with metropolis of the link variables
// retrn 0 if the trial state is rejected and 1 otherwise
int metropolis_for_link(Conf *GC,
                        Geometry const * const geo,
                        GParam const * const param,
                        long r,
                        int i)
  {
  double old_energy, new_energy;
  int old_lambda, new_lambda, pstaple;
  int acc=0;

  #ifdef DEBUG
  long diff_link, diff_plaq;
  double deltaE, deltaEtest;
  #endif

  #ifdef DEBUG
  diff_link=-link(GC, geo, param);
  GC->lambda[r][i]*=-1;
  diff_link+=link(GC, geo, param); //new-old
  GC->lambda[r][i]*=-1;

  diff_plaq=-plaquette(GC, geo, param);
  GC->lambda[r][i]*=-1;
  diff_plaq+=plaquette(GC, geo, param); //new-old
  GC->lambda[r][i]*=-1;

  deltaE=-(param->d_J)*(double)diff_link-(param->d_K)*(double)diff_plaq;
  #endif

  old_lambda=GC->lambda[r][i];

  pstaple=plaqstaples_for_link(GC, geo, r, i);

  old_energy=-(param->d_J)*(double)(old_lambda* GC->phi[r] * GC->phi[nnp(geo, r, i)] );
  old_energy-=param->d_K*(double) (old_lambda*pstaple );

  new_lambda = -old_lambda;

  new_energy=-old_energy;

  #ifdef DEBUG
  deltaEtest=new_energy-old_energy-deltaE;
  if(fabs(deltaEtest)>MIN_VALUE)
    {
    fprintf(stderr, "Error in Metropolis step [deltaEtest=%lg] (%s, %d)\n", deltaEtest, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  if(old_energy>new_energy)
    {
    GC->lambda[r][i] = new_lambda;
    acc=1;
    }
  else if(casuale()< exp(old_energy-new_energy) )
         {
         GC->lambda[r][i] = new_lambda;
         acc=1;
         }

  return acc;
  }


// staples for phi
int staples_for_phi(Conf *GC,
                    Geometry const * const geo,
                    long r)
  {
  int j;
  int ris;

  ris=0;

  for(j=0; j<STDIM; j++)
     {
     ris+=GC->lambda[r][j]*GC->phi[nnp(geo, r, j)];
     ris+=GC->lambda[nnm(geo,r,j)][j]*GC->phi[nnm(geo, r, j)];
     }

    return ris;
    }


// perform an update with metropolis of the phi variable
// retrn 0 if the trial state is rejected and 1 otherwise
int metropolis_for_phi(Conf *GC,
                      Geometry const * const geo,
                      GParam const * const param,
                      long r)
  {
  double old_energy, new_energy;
  int old_phi, new_phi, staple;
  int acc=0;

  #ifdef DEBUG
  long diff_link;
  double deltaE, deltaEtest;
  #endif

  #ifdef DEBUG
  diff_link=-link(GC, geo, param);
  GC->phi[r]*=-1;
  diff_link+=link(GC, geo, param); //new-old
  GC->phi[r]*=-1;

  deltaE=-(param->d_J)*(double)diff_link;
  #endif

  old_phi=GC->phi[r];

  staple=staples_for_phi(GC, geo, r);

  old_energy=-(param->d_J)*old_phi*staple;

  new_phi = -old_phi;

  new_energy=-old_energy;

  #ifdef DEBUG
  deltaEtest=new_energy-old_energy-deltaE;
  if(fabs(deltaEtest)>MIN_VALUE)
    {
    fprintf(stderr, "Error in Metropolis step [deltaEtest=%lg] (%s, %d)\n", deltaEtest, __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }
  #endif

  if(old_energy>new_energy)
    {
    GC->phi[r] = new_phi;
    acc=1;
    }
  else if(casuale()< exp(old_energy-new_energy) )
         {
         GC->phi[r] = new_phi;
         acc=1;
         }

  return acc;
  }


// perform a complete update
void update(Conf * GC,
            Geometry const * const geo,
            GParam const * const param,
            double *acc_link,
            double *acc_site)
   {
   long r, asum_link, asum_site;
   int dir;

   #ifdef OPEN_BC
   long counter=0;
   #endif

   // metropolis on links
   asum_link=0;
   for(r=0; r<param->d_volume; r++)
      {
      for(dir=0; dir<STDIM; dir++)
         {
         #ifndef OPEN_BC
           asum_link+=metropolis_for_link(GC, geo, param, r, dir);
         #else
           if(bcsitep(geo,r, dir)==1)
             {
             asum_link+=metropolis_for_link(GC, geo, param, r, dir);
             counter++;
             }
         #endif
         }
      }
   #ifndef OPEN_BC
     *acc_link=((double)asum_link)*param->d_inv_vol;
     *acc_link/=(double)STDIM;
   #else
     *acc_link=((double)asum_link)/(double)counter;
     *acc_link/=(double)STDIM;
   #endif

   // metropolis on sites
   asum_site=0;
   for(r=0; r<param->d_volume; r++)
      {
      asum_site+=metropolis_for_phi(GC, geo, param, r);
      }
   *acc_site=((double)asum_site)*param->d_inv_vol;

   GC->update_index++;
   }


// staples for gauge
int staples_for_gauge(Conf *GC,
                      Geometry const * const geo,
                      long r)
  {
  int j;
  int ris;

  ris=0;

  for(j=0; j<STDIM; j++)
     {
     ris+=GC->lambda[r][j]*GC->gauge[nnp(geo, r, j)];
     ris+=GC->lambda[nnm(geo,r,j)][j]*GC->gauge[nnm(geo, r, j)];
     }

    return ris;
    }


// perform an update with metropolis of the gauge variable
// retrn 0 if the trial state is rejected and 1 otherwise
int metropolis_for_gauge(Conf *GC,
                         Geometry const * const geo,
                         GParam const * const param,
                         long r)
  {
  double old_energy, new_energy;
  int old_gauge, new_gauge, staple;
  int acc=0;


  staple=staples_for_gauge(GC, geo, r);

  old_gauge=GC->gauge[r];
  old_energy=-(param->d_quench_gamma)*old_gauge*staple;

  new_gauge = -old_gauge;
  new_energy=-old_energy;

  if(old_energy>new_energy)
    {
    GC->gauge[r] = new_gauge;
    acc=1;
    }
  else if(casuale()< exp(old_energy-new_energy) )
         {
         GC->gauge[r] = new_gauge;
         acc=1;
         }

  return acc;
  }


// perform the gauge transformation update and perform measures
// return acceptance rate of update
double glass_evolution_and_meas(Conf *GC,
                                Conf *GC2,
                                GParam const * const param,
                                Geometry const * const geo,
                                FILE *datafilep)
  {
  int err;
  long count, count2, r, nummeas, acc_loc;
  double acc, buffer[4], buffer2[2], *data;

  nummeas=(param->d_quench_sample-param->d_quench_thermal)/param->d_quench_measevery;

  err=posix_memalign((void**) &(data), (size_t) DOUBLE_ALIGN, (size_t) (6*nummeas) * sizeof(double));
  if(err!=0)
    {
    fprintf(stderr, "Problems in allocating the vector for measures! (%s, %d)\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
    }

  acc=0.0;

  // intialize GC2 using GC
  equal_conf(GC2, GC, param);

  count2=0;
  //update both GC and GC2 (gauge variables only)
  for(count=1; count<=param->d_quench_sample; count++)
     {
     acc_loc=0;
     for(r=0; r<param->d_volume; r++)
        {
        acc_loc+=metropolis_for_gauge(GC, geo, param, r);
        acc_loc+=metropolis_for_gauge(GC2, geo, param, r);
        }
     acc+=((double)acc_loc)*param->d_inv_vol*0.5;

     if(count % param->d_quench_measevery ==0 && count > param->d_quench_thermal)
       {
       gauge_apply(GC, geo, param);
       perform_vec_measures_buffer(GC, param, geo, buffer);
       gauge_apply(GC, geo, param);

       perform_overlap_measures_buffer(GC, GC2, param, geo, buffer2);

       data[6*count2+0]=buffer[0]; // vector p=0
       data[6*count2+1]=buffer[1]; // vector pmin
       data[6*count2+2]=buffer[2]; // gauge p=0
       data[6*count2+3]=buffer[3]; // gauge pmin
       data[6*count2+4]=buffer2[0]; // overlap p=0
       data[6*count2+5]=buffer2[1]; // overlap pmin

       count2++;
       }
     }

  acc/=(double)param->d_quench_sample;

  for(r=0; r<6*nummeas; r++)
     {
     fprintf(datafilep, "%.12lf ", data[r]);
     }
  fprintf(datafilep, "\n");
  fflush(datafilep);

  free(data);

  return acc;
  }


#endif
