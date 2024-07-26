INC := -I include/
SRC_OBJS := obj/cs_matrix.o obj/matrix_mul.o \
		obj/csc_io.o obj/radixsort.o obj/transpose.o
TEST_OBJS := obj/matrix_mul_tests.o obj/csc_io_tests.o obj/tests.o \
			 obj/transpose_tests.o

CC = gcc
CFLAGS += -Wall -Wextra -Wpedantic $(INC) -c

all: CFLAGS += -O2 
all: matrixMul

testperf: CFLAGS += -O2
testperf: LDFLAGS += -g
testperf: test

test: CFLAGS += -g
test: $(TEST_OBJS) $(SRC_OBJS)
	$(CC) $(LDFLAGS) $(INC) $^ -o $@

matrixMul: obj/main.o $(SRC_OBJS)
	$(CC) $(LDFLAGS) $(INC) $^ -o $@

obj/main.o: src/main.c include/csc_io.h include/matrix_mul.h include/cs_matrix.h include/transpose.h
	@mkdir -p obj/
	$(CC) $(CFLAGS) $< -o $@

obj/tests.o: tests/tests.c include/matrix_mul_tests.h include/csc_io_tests.h \
				include/transpose_tests.h include/cs_matrix.h include/matrix_mul.h
	@mkdir -p obj/
	$(CC) $(CFLAGS) $< -o $@

obj/cs_matrix.o: src/cs_matrix.c include/cs_matrix.h
	@mkdir -p obj/
	$(CC) $(CFLAGS) $< -o $@

obj/transpose.o: src/transpose.c include/cs_matrix.h include/radixsort.h include/transpose.h
	@mkdir -p obj/
	$(CC) $(CFLAGS) $< -o $@

obj/transpose_tests.o: tests/transpose_tests.c include/transpose_tests.h include/transpose.h include/cs_matrix.h include/radixsort.h include/csc_io.h
	@mkdir -p obj/
	$(CC) $(CFLAGS) $< -o $@

obj/%_tests.o: tests/%_tests.c include/%_tests.h include/%.h include/cs_matrix.h
	@mkdir -p obj/
	$(CC) $(CFLAGS) $< -o $@

obj/%.o: src/%.c include/%.h include/cs_matrix.h
	@mkdir -p obj/
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f main test 
	rm -rf obj/


.PHONY: all clean testperf

