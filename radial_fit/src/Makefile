# Makefile

FILENAME1=fit_Tx_sl_1styear
FILENAME2=fit_sx_vikhlinin_mos1+mos2+pn
FILENAME3=fit_sx_beta_mos1+mos2+pn
FILENAME4=fit_sx_2beta_mos1+mos2+pn

OBJEXT = o
#CXX = g++
#CXXFLAGS = -fopenmp
CXX = icpc
CXXFLAGS = -parallel
CXXLINK = $(CXX) -g -O2 $(CXXFLAGS)  \
	-o $@
normal_OBJECTS1 = ${FILENAME1}.cxx
normal_OBJECTS2 = ${FILENAME2}.cxx
normal_OBJECTS3 = ${FILENAME3}.cxx
normal_OBJECTS4 = ${FILENAME4}.cxx

output1 = ../bin/${FILENAME1}
output2 = ../bin/${FILENAME2}
output3 = ../bin/${FILENAME3}
output4 = ../bin/${FILENAME4}

#INCS = -I${HEADAS}/include \
#	-I/cluster491/software/heasoft/heasoft-6.19/Xspec/src/XSFit/FitMethod/Minuit/minuit2/inc/
#LIBS = -L${HEADAS}/lib  -lcfitsio.5.3.39 -lMinuit2 \
#        -L/usr/lib64 -lgsl

INCS = -I${HEADAS}/include  \
	-I/usr/local/include/Minuit2
LIBS = -L${HEADAS}/lib  -lcfitsio -lMinuit2 \
        -I/usr/local/include -L/usr/lib -lgsl -lgslcblas

all : ${output1} ${output2} ${output3} ${output4}

${output1}: ${FILENAME1}.cxx
	rm -f $@
	$(CXXLINK) $(normal_OBJECTS1) $(INCS) $(LIBS)

${output2}: ${FILENAME2}.cxx
	rm -f $@
	$(CXXLINK) $(normal_OBJECTS2) $(INCS) $(LIBS)

${output3}: ${FILENAME3}.cxx
	rm -f $@
	$(CXXLINK) $(normal_OBJECTS3) $(INCS) $(LIBS)

${output4}: ${FILENAME4}.cxx
	rm -f $@
	$(CXXLINK) $(normal_OBJECTS4) $(INCS) $(LIBS)


clean:
	rm -f ${output1} 
	rm -f ${output2} 
	rm -f ${output3} 
	rm -f ${output4} 
