OMP_EN_VALUE?= 1
OMP_THREAD_NUM ?= 4
export OMP_EN_VALUE
export OMP_THREAD_NUM


all:  
	cd cpp && make && cd ..

.PHONY: test demo

ifeq ($(OMP_EN_VALUE),1)
test: Makefile
	cd cpp && make test && cd ..
	@echo "Using $(OMP_THREAD_NUM) threads"

demo: Makefile
	cd cpp && make demo && cd ..
	@echo "Using $(OMP_THREAD_NUM) threads"
else
test: Makefile
	cd cpp && make test && cd ..

demo: Makefile
	cd cpp && make demo && cd ..
endif

.PHONY: clean
clean:
	cd cpp && make clean && cd ..
	rm -rf 
