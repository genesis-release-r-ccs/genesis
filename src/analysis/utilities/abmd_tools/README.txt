ABMD reweighting tools
======================

This directory contains utility scripts for reweighting GENESIS ABMD simulations.

The basic workflow consists of two steps:

1. Calculate local (unnormalized) reweighting factors with calc_reweight.py.
2. Normalize the weights with calc_normalize_weights.py.

For multiple-walker simulations, calc_reweight.py should first be run
independently for each walker. The resulting log_weight values are then
normalized together over all walkers by calc_normalize_weights.py.


1. Calculate local weights
--------------------------

calc_reweight.py calculates the Tiwary/WT reweighting factor from the
time-dependent ABMD bias.

The script uses two meta files:

  --bias-meta-file
      List of ABMD bias-history (.bp) files.

  --log-meta-file
      List of GENESIS output log or remlog files.

The files in each meta file must be listed in chronological order.


Bias meta file
~~~~~~~~~~~~~~

Example: meta_bias.out

    run_abmd1.bp
    run_abmd2.bp
    run_abmd3.bp
    run_abmd4.bp


Log meta file
~~~~~~~~~~~~~

For a conventional ABMD simulation, meta_log.out may contain normal
GENESIS output files:

    run_abmd1.out
    run_abmd2.out
    run_abmd3.out
    run_abmd4.out

For MW-ABMD or BE-ABMD, the same type of meta file may instead contain
remlog files:

    run_remlog1.log
    run_remlog2.log
    run_remlog3.log
    run_remlog4.log

The log files are interpreted from their contents. Both normal GENESIS
output files and remlog files are accepted as long as they contain
INFO: records including STEP, TIME, and ABMD_BIAS.

The INFO: record at STEP = 0 is ignored.


Example
~~~~~~~

    python calc_reweight.py \
        --bias-meta-file meta_bias.out \
        --log-meta-file meta_log.out \
        -t 300 \
        -wt 3000 \
        -o weight_w1.dat

Options used above:

  -t
      Simulation temperature in K.

  -wt
      GENESIS well-tempered temperature DeltaT in K.

  -o
      Output weight file.

The default output contains:

    # snapshot log_weight
    1 ...
    2 ...
    3 ...

The log_weight values are not normalized at this stage.


2. Repeat for each walker
-------------------------

For MW-ABMD, calculate the local weights independently for all walkers.

For example:

    python calc_reweight.py \
        --bias-meta-file meta_bias.out \
        --log-meta-file meta_log_w1.out \
        -t 300 -wt 3000 \
        -o weight_w1.dat

    python calc_reweight.py \
        --bias-meta-file meta_bias.out \
        --log-meta-file meta_log_w2.out \
        -t 300 -wt 3000 \
        -o weight_w2.dat

Repeat this for all walkers.

Because MW-ABMD walkers share the same ABMD bias history,
--bias-meta-file can be common to all walkers, whereas
--log-meta-file is prepared separately for each walker.


3. Prepare a meta file for the local weights
--------------------------------------------

Create a meta file listing the local-weight files.

Example: meta_weight.out

    weight_w1.dat
    weight_w2.dat
    weight_w3.dat
    weight_w4.dat

For a conventional single-walker ABMD simulation, the meta file contains
only one file:

    weight.dat


4. Normalize the weights
------------------------

Run:

    python calc_normalize_weights.py meta_weight.out

By default, the normalized files are written to:

    global_weights/

For example:

    global_weights/
      weight_w1_global.dat
      weight_w2_global.dat
      weight_w3_global.dat
      weight_w4_global.dat
      global_weight_summary.dat

calc_normalize_weights.py reads all log_weight values from all files listed
in meta_weight.out and computes one global normalization constant over all
walkers and snapshots.

For every walker k and snapshot t:

    global_weight(k,t)
      = exp(log_weight(k,t) - log_Z_global)

where

    log_Z_global
      = log sum_{k,t} exp(log_weight(k,t))

Therefore:

    sum_{k,t} global_weight(k,t) = 1

For a multiple-walker simulation, the weights within an individual walker
do not generally sum to 1. The normalization is performed over the complete
set of walkers.

For a conventional single-walker ABMD simulation, the same procedure can be
used with a meta file containing only one local-weight file.
