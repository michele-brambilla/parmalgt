#!/bin/sh

Lato=$1
Tau=$2
Damocle=$3
Conf="wil_nf0_"$Lato"x"$Lato"x"$Lato"x"$Lato"_a0.05_tg"$Tau"_o6."
Plaq="plaq_wil_nf0_"$Lato"x"$Lato"x"$Lato"x"$Lato"_a0.05_tg"$Tau"_o6."
Seed=$RANDOM

echo $Conf"dat"
echo $Plaq"txt"

#./Quenched.exe -l $Lato -g $Tau -d $Damocle -r $Seed -t 0
#cp result/$Conf"dat" result/$Conf"0"
#mv result/plaquette/$Plaq"txt" result/plaquette/$Plaq"0"

for ((i=19;i<=30;i++))
do
    Seed=$RANDOM
    ./Quenched.exe -l $Lato -g $Tau -d $Damocle -r $Seed -t 1
    cp result/$Conf"dat" result/$Conf$i
    mv result/plaquette/$Plaq"txt" result/plaquette/$Plaq$i
done
