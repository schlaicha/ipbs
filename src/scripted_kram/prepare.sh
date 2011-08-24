. params.sh

if [ -e execute.sh ]; then
  rm execute.sh
fi

# for d in $dist; do
for x in {305..600..5}; do
  d="$(echo "scale=1; ${x}. / 10." | bc)"
  dir=two_colloid_dist_${d}
  if [ -e $dir ]; then
    echo "$dir already exists"
  else
    mkdir $dir
  fi
  cd $dir
  sed -e "s/@distance@/$d/g" ../mesh.geo_template > mesh.geo
  # sed -e "s/@distance@/$d/g" ../mesh.geo_template > .tmp_mesh.geo
  # if [ "$(echo ${d} / .5 | bc)"  -lt 1 ]; then
  #   sed -e "s/@lc_center@/0.1/g" .tmp_mesh.geo > mesh.geo
  # elif [ "$(echo ${d} / 2. | bc)"  -lt 1 ]; then
  #   sed -e "s/@lc_center@/0.25/g" .tmp_mesh.geo > mesh.geo
  # elif [ "$(echo ${d} / 5. | bc)"  -lt 1 ]; then
  #   sed -e "s/@lc_center@/0.5./g" .tmp_mesh.geo > mesh.geo
  # elif [ "$(echo ${d} / 10. | bc)"  -lt 1 ]; then
  #   sed -e "s/@lc_center@/1.0/g" .tmp_mesh.geo > mesh.geo
  # else
  #   sed -e "s/@lc_center@/2.0/g" .tmp_mesh.geo > mesh.geo
  # fi
  # rm .tmp_mesh.geo
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
done
