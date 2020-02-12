#!/bin/bash
name=`basename ${1} .bed`

macs2 callpeak -t ${1} --shift -25 --extsize 50 -f BED --name ${name} --nomodel --keep-dup all -g hs
