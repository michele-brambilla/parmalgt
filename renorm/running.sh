#!/bin/bash

echo "Momenti completati:" 
grep -v ^$ running.log | grep -v [a-z] | wc -l
echo ""
echo "su un totale di:"
cat momREF.txt | wc -l
echo ""