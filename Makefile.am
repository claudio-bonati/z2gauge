bin_PROGRAMS = z2broken

conf = ./include/conf.h \
       ./lib/conf_def.c \
       ./lib/conf_meas.c \
       ./lib/conf_upd.c 

dSFMT = ./include/dSFMT-common.h \
        ./include/dSFMT.h \
        ./include/dSFMT-params11213.h \
        ./include/dSFMT-params1279.h \
        ./include/dSFMT-params132049.h \
        ./include/dSFMT-params19937.h \
        ./include/dSFMT-params216091.h \
        ./include/dSFMT-params2203.h \
        ./include/dSFMT-params4253.h \
        ./include/dSFMT-params44497.h \
        ./include/dSFMT-params521.h \
        ./include/dSFMT-params86243.h \
        ./include/dSFMT-params.h \
        ./lib/dSFMT.c

endian = ./include/endianness.h  \
         ./lib/endianness.c 

geometry = ./include/geometry.h \
           ./lib/geometry.c

gparam = ./include/gparam.h \
         ./lib/gparam.c

random = ./include/myrngparam.h \
         ./include/random.h \
         ./lib/random.c 

z2broken_SOURCES=./src/z2broken.c \
                 $(conf) \
                 $(dSFMT) \
                 $(endian) \
                 $(geometry) \
                 $(gparam) \
                 ./include/macro.h \
                 $(random) 

AM_CFLAGS = -O2 -Wall -Wextra -Werror -pedantic -Wconversion

CLEANFILES = dati.dat log.dat dati.log conf.dat* template_input.in 
