#!/bin/sh

Lato=$1
Tau=$2
Damocle=$3

Ord=`grep allocORD choices.h | awk '{print $3}'`

Ntt=`grep ntt choices.h | head -1 | awk '{print $3}'`
Ntx=`grep ntx choices.h | head -1 | awk '{print $3}'`
Nty=`grep nty choices.h | head -1 | awk '{print $3}'`
Ntz=`grep ntz choices.h | head -1 | awk '{print $3}'`
Nthr=$((Ntt*Ntx*Nty*Ntz))

echo Nthr = $Nthr

Conf="iwa_nf4_"$Lato"x"$Lato"x"$Lato"x"$Lato"_a0.05_tg"$Tau"_tf"$Tau"_o"$Ord"."
Plaq="plaq_iwa_nf4_"$Lato"x"$Lato"x"$Lato"x"$Lato"_a0.05_tg"$Tau"_tf"$Tau"_o"$Ord"."

echo $Conf"dat"
echo $Plaq"txt"


cat Wil.cfg.templ > Wil.cfg

for i in $(seq 1 $Nthr);
  do
  echo seed $RANDOM$RANDOM$RANDOM >> Wil.cfg
done

./Unquenched.exe -l $Lato -g $Tau -f $Tau -d $Damocle -t 0 -n 2

for o in $(seq 4 2 $Ord);
  do

  cat Wil.cfg.templ > Wil.cfg

  for i in $(seq 1 $Nthr);
    do
    echo seed $RANDOM$RANDOM$RANDOM >> Wil.cfg
  done
  
  ./Unquenched.exe -l $Lato -g $Tau -f $Tau -d $Damocle -t 1 -n $o

done


for ((i=100;i<=110;i++))
do
  
  cat Wil.cfg.templ > Wil.cfg
  
  for i in $(seq 1 $Nthr);
    do
    echo seed $RANDOM$RANDOM$RANDOM >> Wil.cfg
  done
  
  ./Unquenched.exe -l $Lato -g $Tau -f $Tau -d $Damocle -t 1 -n $Ord
  cp result/$Conf"dat" result/$Conf$i
  mv result/plaquette/$Plaq"txt" result/plaquette/$Plaq$i
  mv $Conf"trU" $Conf"trU."$i
  mv $Conf"nor" $Conf"nor."$i
  mv $Conf"log" $Conf"log."$i
  ./tools/./NewWilLoop.exe -l $Lato -c ./result/$Conf$i -o ./result/nloop_$Conf$i >& outl0
done
