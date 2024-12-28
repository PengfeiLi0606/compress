#!/bin/sh
ln -s ../cmplib/aarrays.f .  > /dev/null 2>&1
ln -s ../cmplib/data_array.f .  > /dev/null 2>&1
ln -s ../cmplib/inc/ .  > /dev/null 2>&1
gfortran -w -ffpe-trap=invalid,zero,overflow aarrays.f compress.f \
  ../../lib/cmplib.a ../../lib/libutils.a -o ../../compress
rm -f *.o *.mod aarrays.f data_array.f inc
