MUMPS_I = -I../MUMPS_4.10.0_scotch/include -I../MUMPS_4.10.0_scotch/libseq -I../MUMPS_4.10.0_scotch/scotch_5.1.12_esmumps/include
PKG_LIBS = -L../MUMPS_4.10.0_scotch/lib -ldmumps -lmumps_common -lpord  -L../MUMPS_4.10.0_scotch/libseq -lmpiseq -lpthread -L../MUMPS_4.10.0_scotch/scotch_5.1.12_esmumps/lib -lesmumps -lscotch -lscotcherr
