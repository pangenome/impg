#!/usr/bin/env bash
cd /home/erikg/impg-genome-inference
source build-env.sh
taskset -c 252-255 cargo build --release --example router_spell_diag > build-spelldiag3.log 2>&1
echo $? > build-spelldiag3.exit
touch build-spelldiag3.done
