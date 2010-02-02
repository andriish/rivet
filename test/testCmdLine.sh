#!/bin/bash

cp testApi.hepmc file2.hepmc
(rivet -a D0_2008_S7554427 testApi.hepmc file2.hepmc | grep -q "20 events") || exit 1
(cat testApi.hepmc | rivet -a D0_2008_S7554427 | grep -q "10 events") || exit 1
mkfifo fifo.hepmc
(cat testApi.hepmc > fifo.hepmc & rivet -a D0_2008_S7554427 fifo.hepmc | grep -q "10 events") || exit 1
rm -f fifo.hepmc
rm -f file2.hepmc
