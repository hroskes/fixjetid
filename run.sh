#!/bin/bash

set -euo pipefail

cd /afs/cern.ch/work/h/hroskes/public/CJLST/CMSSW_10_2_15
eval $(scram ru -sh)
cd /eos/user/h/hroskes/CJLST/190821/fixjetid
mkdir -p $(dirname /eos/home-h/hroskes/CJLST/190821_fixjetid/$1)
./fixjetid.py /eos/home-h/hroskes/CJLST/190821/$1 /eos/home-h/hroskes/CJLST/190821_fixjetid/$1
