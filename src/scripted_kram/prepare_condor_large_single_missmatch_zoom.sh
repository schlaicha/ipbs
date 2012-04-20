#!/bin/bash

date=$(date +"%F_%T")
rundir=run_$date
runpath=/data/data4/schlaich/data/iPBS/$rundir
execpath=/data/data4/schlaich/progs/dune-2.1_seq/ipbs/src
scriptpath=$execpath/scripted_kram
returnpath=$(pwd)

echo $runpath

if [ -e $runpath ]; then
  echo "Runpath existiert bereits!!"
else
  mkdir $runpath
fi

if [ -e execute.sh ]; then
  rm execute.sh
fi

cd $runpath
for x in `seq 1 25`; do
  d="$(echo "scale=3; 1.02^${x}.*132.0" | bc)"
  echo $d
  dir=two_colloid_dist_${d}
  if [ -e $dir ]; then
    echo "$dir already exists"
  else
    mkdir $dir
  fi
  cd $dir
  sed -e "s/@distance@/$d/g" $scriptpath/2sphere_sym_large.geo_template > mesh.geo
  sed -e "s|@path@|$runpath/$dir|g" $scriptpath/config_large.cfg_template > config.cfg
 
  
  cd ..
  echo "cd $runpath/$dir" >> execute.sh
  echo "condor_run \"gmsh -2 mesh.geo; $execpath/ipbs config.cfg\" 2>err.txt 1>out.txt < /dev/null &" >> $runpath/execute.sh
  
  echo "cd .." >> execute.sh
  echo "head -n 1 $dir/forces.dat | awk '{print \"$d \" \$2 \" \" \$3}' >> force.dat" >> analyze.sh
done
cd $runpath
. execute.sh
echo "Jobs successfully started."

cd $returnpath
