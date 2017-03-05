#!/bin/bash
g95 -c -I ../../ logrid.f90
g95 -c -I ../../ my_atomic.f90
g95 -c -I ../../ ps_hgh.f90
#g95 -c -I ../../ test_hgh.f90
g95 -c -I ../../ hgh_info.f90
g95 -c -I ../../ dump_logrid.f90
