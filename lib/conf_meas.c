#ifndef CONF_MEAS_C
#define CONF_MEAS_C

#include"../include/macro.h"

#include<complex.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

#include"../include/gparam.h"
#include"../include/geometry.h"
#include"../include/conf.h"

// computation of the plaquette in position r and positive directions i,j
int plaquette_single(Conf const * const GC,
                     Geometry const * const geo,
                     long r,
                     int i,
                     int j)
   {

//
//       ^ i
//       |  (3)
//       +---<---+
//       |       |
//   (4) V       ^ (2)
//       |       |
//       +--->---+---> j
//       r  (1)
//

   int ris;

   ris = GC->lambda[r][j];  // (1)
   ris *= GC->lambda[nnp(geo, r, j)][i]; // (2)
   ris *= GC->lambda[nnp(geo, r, i)][j]; // (3)
   ris *= GC->lambda[r][i];

   return ris;
   }


long plaquette(Conf const * const GC,
                 Geometry const * const geo,
                 GParam const * const param)
   {
   long r;
   long ris=0;

   for(r=0; r<(param->d_volume); r++)
      {
      long tmp;
      int i, j;

      i=0;
      tmp=0;
     
      for(i=0; i<STDIM; i++)
         {
         for(j=i+1; j<STDIM; j++)
            {
            tmp+= plaquette_single(GC, geo, r, i, j);
            }
         }

      ris+=tmp;
      }

   return ris;
   }


// compute the average value of phi_x lambda_{x,mu} phi_{x+mu}
long link(Conf const * const GC,
          Geometry const * const geo,
          GParam const * const param)
  {
  int i;
  long r;
  long ris=0;

  for(r=0; r<(param->d_volume); r++)
     {
     for(i=0; i<STDIM; i++)
        {
        ris+= (GC->phi[r])*(GC->lambda[r][i])*(GC->phi[nnp(geo, r, i)]);
        }
     }

  return ris;
  }


void perform_measures(Conf *GC,
                      GParam const * const param,
                      Geometry const * const geo,
                      FILE *datafilep)
   {
   long plaq, avlink;

   avlink=link(GC, geo, param);
   plaq=plaquette(GC, geo, param);

   fprintf(datafilep, "%ld %ld ", avlink, plaq);
   fprintf(datafilep, "\n");

   fflush(datafilep);
   }


// apply the gauge transformation to phi and lambda variables
void gauge_apply(Conf *GC,
                 Geometry const * const geo,
                 GParam const * const param)
  {
  int i;
  long r;

  for(r=0; r<param->d_volume; r++)
     {
     (GC->phi[r])*=(GC->gauge[r]);

     for(i=0; i<STDIM; i++)
        {
        GC->lambda[r][i]*=(GC->gauge[r]*GC->gauge[nnp(geo, r, i)]);
        }
     }
  }


// perform vector-related measures and save results in a buffer
void perform_vec_measures_buffer(Conf *GC,
                                 GParam const * const param,
                                 Geometry const * const geo,
                                 double buffer[2])
   {
   (void) geo; // kept just for consistency with other measures

   int coord[STDIM];
   long r;
   const double p = 2.0*PI/(double)param->d_size[1];
   double V;
   double complex Vp;

   // V =sum_x phi_x
   // Vp=sum_x e^{ipx}phi_x
   V=0.0;
   Vp=0.0+0.0*I;

   for(r=0; r<(param->d_volume); r++)
      {
      V+=(double)GC->phi[r];

      si_to_cart(coord, r, param);

      Vp+=((double complex)GC->phi[r]) * cexp(I*((double)coord[1])*p);
      }

   buffer[0]=V*V*param->d_inv_vol; // tildeG0
   buffer[1]=cabs(Vp)*cabs(Vp)*param->d_inv_vol; // tildeGminp
   }


// performe overlap measures
void perform_overlap_measures_buffer(Conf const * const GC,
                                     Conf const * const GC2,
                                     GParam const * const param,
                                     Geometry const * const geo,
                                     double* buffer)
  {
  long int r, overlap;

  (void) geo; // geo is used only for consistency with other cases
              // so it is possibile to add \xi computation without changes

  overlap=0.0;
  for(r=0; r<param->d_volume; r++)
     {
     overlap+=GC->gauge[r]*GC2->gauge[r];
     }

  *buffer=(double)overlap*(double)overlap*param->d_inv_vol;
  }


#endif
