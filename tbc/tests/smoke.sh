#!/bin/bash
# Smoke test: the CLI must run for every interface class without crashing.
set -e
cd "$(dirname "$0")/.."
run() { echo "== $*"; python tbc.py "$@" --T 300 2>&1 | grep -v -i warn | tail -1; }
run AlN/0001 GaN/0001
run GaN/0001 SiC4H/0001
run AlN/0001 SiC4H/0001
run Si/100 Ge/100
run Si/110 Si/110
run GaAs/100 AlAs/100
run Al/111 Al2O3/0001 --metal Al
run Cu/111 C/111 --metal Cu --polycrystalline
run Au PMMA --metal Au --Kn 3.9e19
run SiO2 Si/100
run GaN/0001 Si/111
run Ti/0001 GaN/0001 --metal Ti
echo "== stack"; python tbc.py --stack "Al/111 | SiO2:1.5 | Si/100" --metal Al --T 300 2>&1 | grep -v -i warn | tail -1
run Pt/111 Si/111 --metal Pt
run Ag PMMA --metal Ag
run Al/111 Cu/111
run GaN/0001 diamond/111
run Si/100 Ge/100 --local-equilibrium
echo "== stack Au/Ti/PMMA"; python tbc.py --stack "Au | Ti:2.4 | PMMA:0 | SiO2:2.05 | Si/100" --T 300 2>&1 | grep -v -i warn | tail -1
echo SMOKE OK
