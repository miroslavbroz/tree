
obj = \
  qsort.o \
  node.o \
  tree.o \

test: test.f90 $(obj)
	gfortran $(obj) -o $@ $<

$(obj): %.o:%.f90
	gfortran -c -o $@ $<

