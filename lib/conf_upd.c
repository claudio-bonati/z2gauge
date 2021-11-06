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
  diff_link=link(GC, param);
  GC->lambda[r][i]*=-1;
  diff_link-=link(GC, param); //old-new
  GC->lambda[r][i]*=-1;

  diff_plaq=plaquette(GC, geo, param);
  GC->lambda[r][i]*=-1;
  diff_plaq-=plaquette(GC, geo, param); //old-new
  GC->lambda[r][i]*=-1;

  deltaE=(param->d_J)*(double)diff_link+(param->d_K)*(double)diff_plaq;
  #endif

  old_lambda=GC->lambda[r][i];

  pstaple=plaqstaples_for_link(GC, geo, r, i);

  old_energy=-(param->d_J)*(double)old_lambda;
  old_energy-=param->d_K*(double)old_lambda*(double)pstaple;

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


// perform a complete update
void update(Conf * GC,
            Geometry const * const geo,
            GParam const * const param,
            double *acc_link)
   {
   long r, asum_link;
   int dir;

   // metropolis on links
   asum_link=0;
   for(r=0; r<param->d_volume; r++)
      {
      for(dir=0; dir<STDIM; dir++)
         {
         asum_link+=metropolis_for_link(GC, geo, param, r, dir);
         }
      }
   *acc_link=((double)asum_link)*param->d_inv_vol;
   *acc_link/=(double)STDIM;

   GC->update_index++;
   }

#endif
