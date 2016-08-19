all : libs tests

libs : psp_lib MatrixSwitch_lib libOMM_lib tomato_lib

tests : psp_test MatrixSwitch_examples libOMM_examples

clean : clean_libs clean_tests

clean_libs : clean_psp clean_MatrixSwitch clean_libOMM clean_tomato

clean_tests : clean_psp_test clean_MatrixSwitch_examples clean_libOMM_examples

psp_lib :
	cd psp; \
	make; \
	make install

MatrixSwitch_lib : psp_lib
	cd MatrixSwitch/src; \
	make; \
	make install

libOMM_lib : MatrixSwitch_lib psp_lib
	cd libOMM/src; \
	make; \
	make install

tomato_lib : MatrixSwitch_lib psp_lib
	cd tomato/src; \
	make; \
	make install

psp_test : psp_lib
	cd psp/test; \
	make

MatrixSwitch_examples : MatrixSwitch_lib
	cd MatrixSwitch/examples; \
	make

libOMM_examples : libOMM_lib
	cd libOMM/examples; \
	make

clean_psp :
	cd psp; \
	make clean

clean_MatrixSwitch :
	cd MatrixSwitch/src; \
	make clean

clean_libOMM :
	cd libOMM/src; \
	make clean

clean_tomato :
	cd tomato/src; \
	make clean

clean_psp_test :
	cd psp/test; \
	make clean

clean_MatrixSwitch_examples :
	cd MatrixSwitch/examples; \
	make clean

clean_libOMM_examples :
	cd libOMM/examples; \
	make clean
