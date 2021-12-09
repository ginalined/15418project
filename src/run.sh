#!/bin/sh
clear;

cd input;
python defineSpheres.py;
python sphere2tri.py;
rm spheres.txt;
# cp spheres.inp ../../../VCOLLIDE/demos/nbody/input/obj20.inp;
cd ../trans; 
python generateTrans.py;
# cp rand.inp ../../../VCOLLIDE/demos/nbody/trans;
cd ../;

./cudaCollide input/spheres.inp trans/y2batch.inp;

# rm input/*.inp;
# rm trans/*.inp;