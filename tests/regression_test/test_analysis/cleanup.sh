#!/bin/sh

rm ./*/*/out
rm ./*/*/err
rm ./*/*/log

rm ./*/*/*/out
rm ./*/*/*/err
rm ./*/*/*/log

rm test_remd_convert/APP/*.log
rm test_remd_convert/APP/*.dcd
rm test_kmeans_clustering/BPTI/output*
rm test_flccrd_analysis/BPTI/output*
rm test_eigmat_analysis/APP/output*
rm test_rpath_generator/AdK/*.rst
rm test_rpath_generator/AdK/*.pdb
rm test_pathcv_analysis/path/*.pathcv
rm test_*/*/*.dat
rm test_spheres_generator/BPTI/*.*
rm test_add_ions/BPTI/*.*
rm test_add_metabolites/BPTI/*.*
rm test_ca_fitting/BPTI/*.*
rm test_collision_canceller/BPTI/*.*
rm test_solvate/BPTI/*.*
rm test_spheres_generator/BPTI/*.*
rm trajectories/BPTI_cssb/BPTI/*.*
