CC=clang++

bin/igrep: obj/kernel.o obj/main.o
	${CC} -o $@ $^ -pthread -L${BOOST_ROOT}/lib -lboost_thread -lboost_system -lboost_filesystem -lboost_iostreams -lboost_date_time -L${CUDA_ROOT}/lib64 -lcudart -L${POCO_ROOT}/lib/Linux/x86_64 -lPocoFoundation -lPocoNet -L${MONGODBCXXDRIVER_ROOT} -lmongoclient

obj/%.o: src/%.cu
	nvcc -o $@ $< -c -gencode arch=compute_35,code=sm_35 -I${CUDA_ROOT}/samples/common/inc #-maxrregcount=N -Xptxas=-v

obj/%.o: src/%.cpp
	${CC} -o $@ $< -c -std=c++11 -O3 -Wall -Wno-deprecated-declarations -Wno-unused-local-typedef -I${BOOST_ROOT} -I${MONGODBCXXDRIVER_ROOT}/src -I${POCO_ROOT}/Foundation/include -I${POCO_ROOT}/Net/include -I${CUDA_ROOT}/include -I${CUDA_ROOT}/samples/common/inc

clean:
	rm -f bin/igrep obj/*.o
