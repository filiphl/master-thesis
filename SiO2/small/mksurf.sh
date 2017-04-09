for f in restartFiles/*
do lmp_mpi -in continue.in -var filename $f
done
