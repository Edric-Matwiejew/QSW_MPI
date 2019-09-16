#python3 -m numpy.f2py -h --overwrite-signature fMPI.pyf -m fMPI fMPI.f90
#python -m numpy.f2py --f90exec="mpifort" --f90flags="-fbounds-check" -c fMPI.pyf ../sparse.f90 ../one_norms.f90 ../expm.f90 ../operators.f90 -m fMPI fMPI.f90
#python3 -m numpy.f2py -h --overwrite-signature foperators.pyf -m foperators foperators.f90
#python -m numpy.f2py --f90exec="mpifort" --f90flags="-fbounds-check" -c ../sparse.f90 ../operators.f90  -m foperators foperators.f90
python3 -m numpy.f2py --f90exec="mpifort" --f90flags="-O3 -mtune=native -ffast-math -funroll-loops" -lgomp -c fMPI.pyf ../sparse.f90 ../one_norms.f90 ../expm.f90 ../operators.f90 -m fMPI fMPI.f90
#python -m numpy.f2py -h --overwrite-signature foperators.pyf -m foperators foperators.f90
#f2py --f90exec="mpifort" --f90flags="-fopenmp -fbounds-check" -lgomp -c ../sparse.f90 ../operators.f90  -m foperators foperators.f90
python3 -m numpy.f2py --f90exec="mpifort" --f90flags="-O3 -mtune=native -ffast-math -funroll-loops" -lgomp -c foperators.pyf ../sparse.f90 ../operators.f90 -m foperators foperators.f90
#
#python3 -m numpy.f2py --f90exec="mpifort" --f90flags="-g" -lgomp -c fMPI.pyf ../sparse.f90 ../one_norms.f90 ../expm.f90 ../operators.f90 -m fMPI fMPI.f90
#python3 -m numpy.f2py --f90exec="mpifort" --f90flags="-g" -lgomp -c foperators.pyf ../sparse.f90 ../operators.f90 -m foperators foperators.f90
#
#python3 -m numpy.f2py --f90exec="mpifort" --f90flags="-g -fbounds-check" -lgomp -c fMPI.pyf ../sparse.f90 ../one_norms.f90 ../expm.f90 ../operators.f90 -m fMPI fMPI.f90
#python3 -m numpy.f2py --f90exec="mpifort" --f90flags="-g -fbounds-check" -lgomp  -c -m fMPI fMPI.f90
#python3 -m numpy.f2py --f90exec="mpifort" --f90flags="-g -fbounds-check" -lgomp -c foperators.pyf ../sparse.f90 ../operators.f90 -m foperators foperators.f90
#
rm ../../freeqsw/*.so
cp *.so ../../freeqsw/
