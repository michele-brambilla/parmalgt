#! /bin/bash

ROOTD="${HOME}/parmalgt-0.11/"
WORKDIR="L4runs"
DATADIR="${HOME}/dataL4"

mkdir -p $WORKDIR
mkdir -p $DATADIR
cd $WORKDIR

for tau in .005 .01 .02 .03
do
  for s in -1 0 1
  do
    NAME="run_s${s}_tau${tau}"
    mkdir -p $NAME
    cd $NAME
    $ROOTD/configure
    make
    sed "s/SPAR/${s}/" $ROOTD/Quench.cfg > Quench.cfg
    sed -i "s/TAU/${tau}/" Quench.cfg
    sed -i "s/TALLOC/$((5-s))/" Quench.cfg
    sed "s/TAU/${tau}/" $ROOTD/damocle.dag > damocle.dag
    ./Quenched &
    cd ..
  done
done

#rsync data to datadir

while [ 1 ]
do
  sleep 2m
  rsync -au --include "+ */" --include "+ *.bindat" --exclude "- *" . $DATADIR
done
