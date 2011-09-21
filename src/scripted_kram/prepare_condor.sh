#!/bin/bish

if [ -e execute.sh ]; then
  rm execute.sh
fi

c=1
for x in {1..240..4}; do
  d="$(echo "scale=4; (${c})/1000+2" | bc)"
  dir=two_colloid_dist_${d}
  if [ -e $dir ]; then
    echo "$dir already exists"
  else
    mkdir $dir
  fi
  cd $dir
  sed -e "s/@distance@/$d/g" ../mesh.geo_template > mesh.geo
 
  cat ../config.cfg_template > config.cfg
  
  cd ..
  echo "cd $dir" >> execute.sh
  echo "condor_run \"gmsh -2 mesh.geo; ../../iPBS_without config.cfg\" 2>err.txt 1>out.txt < /dev/null &" >> execute.sh
  
  echo "cd .." >> execute.sh
  echo "head -n 1 $dir/forces.dat | awk '{print \"$d \" \$2 \" \" \$3}' >> force.dat" >> analyze.sh
  echo "tail -n 1 $dir/forces.dat | awk '{print \"$d \" \$2 \" \" \$3}' >> force2.dat" >> analyze.sh
  echo "head -n 1 $dir/forces2.dat | awk '{print \"$d \" \$2 \" \" \$3}' >> force3.dat" >> analyze.sh
  echo "tail -n 1 $dir/forces2.dat | awk '{print \"$d \" \$2 \" \" \$3}' >> force4.dat" >> analyze.sh
  c=$(( $c + $x ))
done
