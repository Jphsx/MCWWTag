#!/bin/sh

root -l -b -q 'myplot2.C()'
root -l -b -q 'myplot2.C("hxT",1500.0)'
root -l -b -q 'myplot2.C("hcosljet",1500.0)'
root -l -b -q 'myplot2.C("hmqcosthl",2400.0)'
root -l -b -q 'myplot2.C("hmqcosthl",2400.0,true,0.48,0.15)'

exit

