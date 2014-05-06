CC=g++ -std=c++11 -O3

all: bin/idock bin/idock2

bin/idock2: obj/scoring_function.o obj/box.o obj/quaternion.o obj/io_service_pool.o obj/safe_counter.o obj/receptor.o obj/ligand.o obj/grid_map_task.o obj/monte_carlo_task.o obj/random_forest_test.o obj/main2.o
	${CC} -o $@ $^ -pthread -L${BOOST_ROOT}/lib -lboost_thread -lboost_program_options -lboost_system -lboost_filesystem -lboost_iostreams -lboost_date_time -L${POCO_ROOT}/lib/Linux/x86_64 -lPocoFoundation -lPocoNet -L${MONGODBCXXDRIVER_ROOT}/build/linux2/64/release/use-system-all -lmongoclient

bin/idock: obj/scoring_function.o obj/box.o obj/quaternion.o obj/io_service_pool.o obj/safe_counter.o obj/receptor.o obj/ligand.o obj/grid_map_task.o obj/monte_carlo_task.o obj/random_forest_test.o obj/main.o
	${CC} -o $@ $^ -pthread -L${BOOST_ROOT}/lib -lboost_thread -lboost_program_options -lboost_system -lboost_filesystem -lboost_iostreams -lboost_date_time -L${POCO_ROOT}/lib/Linux/x86_64 -lPocoFoundation -lPocoNet -L${MONGODBCXXDRIVER_ROOT}/build/linux2/64/release/use-system-all -lmongoclient

obj/main2.o: src/main2.cpp
	${CC} -o $@ $< -c -I${BOOST_ROOT} -I${MONGODBCXXDRIVER_ROOT}/src -I${POCO_ROOT}/Foundation/include -I${POCO_ROOT}/Net/include -Wno-deprecated-declarations -Wno-deprecated-register

obj/main.o: src/main.cpp
	${CC} -o $@ $< -c -I${BOOST_ROOT} -I${MONGODBCXXDRIVER_ROOT}/src -I${POCO_ROOT}/Foundation/include -I${POCO_ROOT}/Net/include -Wno-deprecated-declarations -Wno-deprecated-register

obj/%.o: src/%.cpp
	${CC} -o $@ $< -c -I${BOOST_ROOT}

clean:
	rm -f bin/idock bin/idock2 obj/*.o
