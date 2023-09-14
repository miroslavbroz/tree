
obj = node.o qsort.o

tree: tree.f90 $(obj)
	gfortran $(obj) -o $@ $<

$(obj): %.o:%.f90
	gfortran -c -o $@ $<

