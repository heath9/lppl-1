#!/bin/csh
make lppl
crm tmp.out
lppl -f gld.plt -p gld.par -S 80 -v > tmp.out
cmp tmp.out gld.out
