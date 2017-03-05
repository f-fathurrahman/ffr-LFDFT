#!/bin/bash
g95 -I ../../ \
   logrid.f90 \
   my_atomic.f90 \
   ps_hgh.f90 \
   test_hgh.f90 \
   hgh_info.f90 \
   dump_logrid.f90
