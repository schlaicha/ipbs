. params.sh

if [ -e execute.sh ]; then
  rm execute.sh
fi

# for d in $dist; do
for x in {20..110..10}; do
  d="$(echo "scale=1; ${x}. / 1." | bc)"
  dir=plane_dist_${d}
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
  echo "cd .." >> execute.sh
  echo "head -n 1 $dir/forces.dat | awk '{print \"$d \" \$2 \" \" \$3}' >> force.dat" >> analyze.sh
  echo "tail -n 1 $dir/forces.dat | awk '{print \"$d \" \$2 \" \" \$3}' >> force2.dat" >> analyze.sh
  echo "head -n 1 $dir/forces2.dat | awk '{print \"$d \" \$2 \" \" \$3}' >> force3.dat" >> analyze.sh
  echo "tail -n 1 $dir/forces2.dat | awk '{print \"$d \" \$2 \" \" \$3}' >> force4.dat" >> analyze.sh
  echo "tail -n 1  $dir/surface_potential.dat | awk '{print \"$d \" \$2}' >> surface_potential.dat" >> analyze.sh
done
