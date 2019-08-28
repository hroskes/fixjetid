#!/bin/bash

set -euo pipefail

cd /work-zfs/lhc/CJLSTtrees/190821/fixjetid/CMSSW_10_2_5
eval $(scram ru -sh)
cd /work-zfs/lhc/CJLSTtrees/190821/fixjetid/

./fixjetid.py "$@"
