#!/bin/sh
ADDR=$(cd "$(dirname "$0")"; pwd)
echo "Working directory: $ADDR"

if [ ! -d ${ADDR}/bin/ ]; then
	mkdir ${ADDR}/bin
fi

#ALC
echo "Checking and compiling ALC"
if [ ! -d ${ADDR}/ApproxLatCount/ ]; then
	unzip ApproxLatCount.zip
	cd ${ADDR}/ApproxLatCount
	make
	cd ../
	mv ${ADDR}/ApproxLatCount/ApproxLatCount ${ADDR}/bin/.
fi

#volce main program
echo "Compiling sharpSMT ..."
cd ${ADDR}
make depend
make
