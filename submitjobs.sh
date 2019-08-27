#!/bin/bash

set -euo pipefail

for a in Data_2016  Data_2017  Data_2018  MC_2016  MC_2016_anomalous  MC_2017  MC_2017_anomalous  MC_2018  MC_2018_anomalous; do
  for b in $(ls /eos/home-h/hroskes/CJLST/190821/$a | grep -vi Chunk); do
    mkdir /eos/home-h/hroskes/CJLST/190821_fixjetid/$a/$b -p
    condor_submit 
