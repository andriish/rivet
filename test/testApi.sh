#! /usr/bin/env bash

rivet-build ../analyses/pluginMC/EXAMPLE.cc ../analyses/pluginMC/MC_JETS.cc ../analyses/pluginMC/EXAMPLE_CUTS.cc ../analyses/pluginMC/EXAMPLE_SMEAR.cc
export RIVET_ANALYSIS_PATH=$RIVET_ANALYSIS_PATH:$PWD

exec ./testApi "$srcdir/testApi.hepmc"
