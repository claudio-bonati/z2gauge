#ifndef CONF_MEAS_C
#define CONF_MEAS_C

#include"../include/macro.h"

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


// compute the average value of Re[ phi_x^{dag} lambda_{x,mu} phi_{x+mu} ]
long link(Conf const * const GC,
          GParam const * const param)
  {
  int i;
  long r;
  long ris=0;

  for(r=0; r<(param->d_volume); r++)
     {
     for(i=0; i<STDIM; i++)
        {
        ris+= (GC->lambda[r][i]);
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

   avlink=link(GC, param);
   plaq=plaquette(GC, geo, param);

   fprintf(datafilep, "%ld %ld ", avlink, plaq);
   fprintf(datafilep, "\n");

   //fflush(datafilep);
   }


#endif
