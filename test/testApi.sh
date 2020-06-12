#! /usr/bin/env bash

bash $top_builddir/bin/rivet-build $top_srcdir/analyses/pluginMC/EXAMPLE.cc $top_srcdir/analyses/pluginMC/MC_JETS.cc $top_srcdir/analyses/pluginMC/EXAMPLE_CUTS.cc $top_srcdir/analyses/pluginMC/EXAMPLE_SMEAR.cc
export RIVET_ANALYSIS_PATH=$RIVET_ANALYSIS_PATH:$PWD

exec ./testApi "$srcdir/testApi.hepmc"
