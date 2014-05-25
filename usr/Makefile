CC=clang++

all: bin/encode bin/usr

bin/encode: src/encode.cpp
	${CC} -o $@ $^ -std=c++11 -O2

bin/usr: src/usr.cpp
	${CC} -o $@ $^ -std=c++11 -O2 -mavx -Wno-deprecated-declarations -Wno-deprecated-register -I${BOOST_ROOT} -I${MONGODBCXXDRIVER_ROOT}/src -I${POCO_ROOT}/Foundation/include -I${POCO_ROOT}/Net/include -L${BOOST_ROOT}/lib -lboost_system -lboost_filesystem -lboost_iostreams -lboost_date_time -L${POCO_ROOT}/lib/Linux/x86_64 -lPocoFoundation -lPocoNet -L${MONGODBCXXDRIVER_ROOT}/build/linux2/64/release/use-system-all -lmongoclient

clean:
	rm -f bin/*
