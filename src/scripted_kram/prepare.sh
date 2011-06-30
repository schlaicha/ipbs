. params.sh

if [ -e execute.sh ]; then
  rm execute.sh
fi

for d in $dist; do
  echo $d
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
  echo "gmsh -2 mesh.geo" >> execute.sh 
  echo "../../iPBS config.cfg" >> execute.sh
  echo "# Forces acting on particle i (z and r component):" > $dir/forces.dat
  echo "cd .." >> execute.sh
  echo "tail -n 2 $dir/forces.dat | awk '{print \"$d \" \$2 \" \" \$3}' | head -n 1 >> force.dat" >> analyze.sh
  echo "tail -n 4 $dir/forces.dat | awk '{print \"$d \" \$2 \" \" \$3}' | tail -n 3 >> force2.dat" >> analyze.sh
done
