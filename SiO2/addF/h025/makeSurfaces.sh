#!/bin/bash

for f in restartFiles/*
do
  ../../../lammps/src/lmp_mpi -in hold.in -var filename $f
done

for f in surfaceFiles/*
do
  mkdir ${f:0:${#f}-5}
  mv $f ${f:0:${#f}-5}/
done
