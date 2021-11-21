#!/bin/bash

rm -r stats

make clean

module rm gcc
module rm OpenCilk/9.0.1

module load gcc
make all

mkdir stats

echo "Sequential" > ./stats/sequential.txt
for file in ./matrices/*; do
	printf "\n"
	echo "Sequential calculation on matrix: "$file >> ./stats/sequential.txt
	./sequential $file >> ./stats/sequential.txt
	echo "" >> ./stats/sequential.txt
done

echo "Openck" > ./stats/openck.txt
for file in ./matrices/*; do
	echo "Openck calculation on matrix: "$file >> ./stats/openck.txt
	for i in 2 4 5 10 15 20; do
		./openck $file $i >> ./stats/openck.txt	
	done
	echo " " >> ./stats/openck.txt
done

echo "Openmp" > ./stats/openmp.txt
for file in ./matrices/*; do
	echo "Openmp calculation on matrix: "$file >> ./stats/openmp.txt
	for i in 2 4 5 10 15 20; do
		./openmp $file $i >> ./stats/openmp.txt	
	done
	echo " " >> ./stats/openmp.txt
done

echo "Pthreads" > ./stats/pthreads.txt
for file in ./matrices/*; do
	echo "Pthreads calculation on matrix: "$file >> ./stats/pthreads.txt
	for i in 2 4 5 10 15 20; do
		./pthreads $file $i >> ./stats/pthreads.txt	
	done
	echo " " >> ./stats/pthreads.txt
done