CC=clang++ -std=c++11 -O3
NVCC=nvcc -gencode arch=compute_13,code=sm_13 -Xptxas=-v -use_fast_math# -maxrregcount=N

bin/igrep: obj/kernel.o obj/main.o
	${CC} -o $@ $^ -pthread -lboost_thread -lboost_system -lboost_filesystem -lboost_iostreams -lboost_date_time -L${CUDA_ROOT}/lib64 -lcudart -lcurand -L${POCO_ROOT}/lib/Linux/x86_64 -lPocoFoundation -lPocoNet -L${MONGODBCXXDRIVER_ROOT}/build/linux2/64/release/use-system-all -lmongoclient

obj/kernel.o: src/kernel.cu
	${NVCC} -o $@ $< -c -I${CUDA_ROOT}/samples/common/inc

obj/main.o: src/main.cpp 
	${CC} -o $@ $< -c -I${MONGODBCXXDRIVER_ROOT}/src -I${POCO_ROOT}/Net/include -I${CUDA_ROOT}/include -I${CUDA_ROOT}/samples/common/inc -Wno-deprecated-declarations

clean:
	rm -f bin/igrep obj/*.o
