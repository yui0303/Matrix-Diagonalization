
all: 
	cd cpp && make && cd ..

.PHONY: test
test:
	cd cpp && make test && cd ..

.PHONY: clean
clean:
	cd cpp && make clean && cd ..
