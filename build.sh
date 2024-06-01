#!/bin/sh
ADDR=$(cd "$(dirname "$0")"; pwd)
echo "Working directory: $ADDR"

if [ ! -d ${ADDR}/bin/ ]; then
	mkdir ${ADDR}/bin
fi

#ALC
echo "Checking and compiling ALC ..."
if [ ! -d ${ADDR}/ApproxLatCount/ ]; then
	unzip ApproxLatCount.zip
	cd ${ADDR}/ApproxLatCount
	make
	cd ../
	mv ${ADDR}/ApproxLatCount/ApproxLatCount ${ADDR}/bin/.
fi

#z3
echo "Checking and installing Z3 ..."
if [ ! -f ${ADDR}/z3-master/lib/libz3.so ]; then
	unzip z3-master.zip
	cd ${ADDR}/z3-master
	python scripts/mk_make.py --prefix=${ADDR}/z3-master
	cd ${ADDR}/z3-master/build
	make
	make install
fi

#vinci
echo "Cheking and compiling Vinci ..."
if [ ! -d ${ADDR}/vinci-1.0.5/ ]; then
	unzip vinci-1.0.5.zip
	cd ${ADDR}/vinci-1.0.5
	make
	cd ../
	mv ${ADDR}/vinci-1.0.5/vinci ${ADDR}/bin/vinci
fi

#volce main program
echo "Compiling sharpSMT ..."
cd ${ADDR}
make depend
make
