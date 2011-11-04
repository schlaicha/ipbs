if [ -e execute.sh ]; then
  rm execute.sh
fi

for x in `seq 5 1 25`; do
  d="$(echo "scale=3; 1.15^${x}. / 3. * 40." | bc)"
  echo $d
  dir=two_colloid_dist_${d}
  if [ -e $dir ]; then
    echo "$dir already exists"
  else
    mkdir $dir
  fi
  cd $dir
  sed -e "s/@distance@/$d/g" ../2sphere_sym.geo_template > mesh.geo
  cat ../config.cfg_template > config.cfg
  cd ..
  echo "cd $dir" >> execute.sh
  echo "mpirun -np 4 gmsh -2 mesh.geo" >> execute.sh
  echo "mpirun -np 4 /home/phate/dune/iPBS/src/iPBS config.cfg" >> execute.sh 
  echo "cd .." >> execute.sh
  echo "head -n 1 $dir/forces.dat | awk '{print \"$d \" \$2 \" \" \$3}' >> force.dat" >> analyze.sh
  echo "tail -n 1 $dir/forces.dat | awk '{print \"$d \" \$2 \" \" \$3}' >> force2.dat" >> analyze.sh
  echo "cat $dir/forces2.dat | awk '{print \"$d \" \$2 \" \" \$3}' >> force3.dat" >> analyze.sh
done
