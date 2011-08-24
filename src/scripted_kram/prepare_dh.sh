. params.sh

if [ -e execute.sh ]; then
  rm execute.sh
fi

for x in {1..50..1}; do
  ka="$(echo "scale=4; 1.2^${x}. /100." | bc)"
  lambda="$(echo "scale=4; 1/${ka}" | bc)"      
  dir=sphere_${ka}
  if [ -e $dir ]; then
    echo "$dir already exists"
  else
    mkdir $dir
  fi
  cd $dir
  #sed -e "s/@radius@/$ka/g" ../sphere.geo_template > mesh.geo
  #sigma="$(echo "scale=4; 1/(4*3.14159*${x}^2)" | bc)"
  #sed -e "s/@sigma@/$sigma/g" ../config.cfg_template > config.cfg
  sed -e "s/@lambda@/$lambda/g" ../config.cfg_template > config.cfg
  cat ../sphere.geo_template > mesh.geo
  cd ..
  echo "cd $dir" >> execute.sh
  echo "gmsh -2 mesh.geo" >> execute.sh 
  echo "../../iPBS config.cfg" >> execute.sh
  echo "cd .." >> execute.sh
  echo "tail -n 1  $dir/surface_potential.dat | awk '{print \"$ka \" \$2}' >> surface_potential.dat" >> analyze.sh
done
