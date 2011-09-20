#!/bin/sh

./quench.sh 2 0.005 damocle.dag.1 >& e1 &
./quench.sh 2 0.010 damocle.dag.2 >& e2 &
./quench.sh 2 0.015 damocle.dag.3 >& e3 &