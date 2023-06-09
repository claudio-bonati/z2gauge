#ifndef CONF_H
#define CONF_H

#include"macro.h"

#include<openssl/md5.h>
#include<stdio.h>

#include"gparam.h"
#include"geometry.h"

typedef struct Conf {
  long update_index;

  int **lambda;    // [volume][STDIM]
  } Conf;

// in conf_def.c
void init_conf(Conf *GC,
               GParam const * const param);
void read_conf(Conf *GC,
               GParam const * const param);
void free_conf(Conf *GC,
               GParam const * const param);
void write_conf_on_file_with_name(Conf const * const GC,
                                  GParam const * const param,
                                  char const * const namefile);
void write_conf_on_file(Conf const * const GC,
                        GParam const * const param);
void write_conf_on_file_back(Conf const * const GC,
                             GParam const * const param);
void compute_md5sum_conf(char *res,        // the lenght is 2*MD5_DIGEST_LENGTH
                         Conf const * const GC,
                         GParam const * const param);


// in conf_update.c
int plaqstaples_for_link(Conf *GC,
                         Geometry const * const geo,
                         long r,
                         int i);
int metropolis_for_link(Conf *GC,
                        Geometry const * const geo,
                        GParam const * const param,
                        long r,
                        int i);
void update(Conf * GC,
            Geometry const * const geo,
            GParam const * const param,
            double *acc_link);

// in conf_meas.c
int plaquette_single(Conf const * const GC,
                     Geometry const * const geo,
                     long r,
                     int i,
                     int j);
long plaquette(Conf const * const GC,
                 Geometry const * const geo,
                 GParam const * const param);
long link(Conf const * const GC,
            GParam const * const param);
void perform_measures(Conf *GC,
                      GParam const * const param,
                      Geometry const * const geo,
                      FILE *datafilep);

#endif
