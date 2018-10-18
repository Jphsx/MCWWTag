#!/bin/sh

root -l -b -q 'myplot.C()'
root -l -b -q 'myplot.C("hs")'

root -l -b -q 'myplotp.C()'
root -l -b -q 'myplotp.C("hs")'

exit

