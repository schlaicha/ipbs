. params.sh

if [ -e execute.sh ]; then
  rm execute.sh
fi

for x in {1..23..1}; do
#for x in {10..32..1}; do
  q="$(echo "scale=5; 2^${x}. /250." | bc)"
  dir=sphere_${q}
  if [ -e $dir ]; then
    echo "$dir already exists"
  else
    mkdir $dir
  fi
  cd $dir
  sed -e "s/@radius@/$ka/g" ../sphere.geo_template > mesh.geo
  #sigma="$(echo "scale=4; 1/(4*3.14159*${x}^2)" | bc)"
  sed -e "s/@sigma@/$q/g" ../config.cfg_template > config.cfg
  #sed -e "s/@lambda@/"0.05"/g" mesh.geo.tmp > mesh.geo
  #cat ../sphere.geo_template > mesh.geo
  #cat ../plane.geo_template > mesh.geo
  cd ..
  echo "cd $dir" >> execute.sh
  echo "gmsh -2 mesh.geo" >> execute.sh 
  echo "../../iPBS config.cfg" >> execute.sh
  echo "cd .." >> execute.sh
  echo "tail -n 1  $dir/surface_potential.dat | awk '{print \"$q \" \$2}' >> surface_potential.dat" >> analyze.sh
done
