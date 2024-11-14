#!/usr/bin/env bash
set -o pipefail
git clone https://github.com/ptrebert/sciddo.git assets/bin/sciddo
cd assets/bin/sciddo || exit
python setup.py install
