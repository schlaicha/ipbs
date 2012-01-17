#!/bin/bash

date=$(date +"%F_%T")
rundir=run_$date
runpath=/data/data4/schlaich/data/iPBS/$rundir
execpath=/data/data4/schlaich/progs/dune-2.1/iPBS/src
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
  d="$(echo "scale=3; 1.02^${x}.*80.0" | bc)"
  echo $d
  dir=two_colloid_dist_${d}
  if [ -e $dir ]; then
    echo "$dir already exists"
  else
    mkdir $dir
  fi
  cd $dir
  sed -e "s/@distance@/$d/g" $scriptpath/2sphere_sym_medium.geo_template > mesh.geo
  sed -e "s|@path@|$runpath/$dir|g" $scriptpath/config_medium.cfg_template > config.cfg
 
  
  cd ..
  echo "cd $runpath/$dir" >> execute.sh
  echo "condor_run \"gmsh -2 mesh.geo; $execpath/iPBS config.cfg\" 2>err.txt 1>out.txt < /dev/null &" >> $runpath/execute.sh
  
  echo "cd .." >> execute.sh
  echo "head -n 1 $dir/forces.dat | awk '{print \"$d \" \$2 \" \" \$3}' >> force.dat" >> analyze.sh
  echo "tail -n 1 $dir/forces.dat | awk '{print \"$d \" \$2 \" \" \$3}' >> force2.dat" >> analyze.sh
  echo "cat $dir/forces2.dat | awk '{print \"$d \" \$2}' >> force3.dat" >> analyze.sh
done
cd $runpath
. execute.sh
echo "Jobs successfully started."

cd $returnpath
