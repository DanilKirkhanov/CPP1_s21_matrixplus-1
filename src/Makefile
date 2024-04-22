CC = g++
FLAGS = -Wall -Werror -Wextra -pedantic -std=c++17
LIB = s21_matrix_oop
FILES=$(wildcard s21_*.cpp)
OBJ=$(patsubst %.cpp,%.o, ${FILES})
GCOV_FLAGS=--coverage
T_FILES= test_me.cpp
UNAME := $(shell uname)

ifeq ($(UNAME), Linux)
	OPEN_CMD = xdg-open
	TEST_LIB = -lgtest -lgtest_main -lm -pthread -lrt -lsubunit
endif
ifeq ($(UNAME), Darwin)
	OPEN_CMD = open
	TEST_LIB = -lcheck -lgtest -lgtest_main
endif

all: clean s21_matrix_oop.a test

s21_matrix_oop.a: $(OBJ)
	@ar rc $@ $(OBJ)
	@ranlib $@

$(OBJ): $(FILES)
	$(CC) $(FLAGS) $(GCOV_FLAGS) -c $(FILES)

test: clean $(LIB).a
	$(CC) $(FLAGS) $(GCOV_FLAGS) -o test test.cpp $(LIB).a $(TEST_LIB)
	./test

gcov_report: test clean
	gcc  --coverage test.cpp s21_matrix_oop.cpp -o gcov_report -lgtest -lstdc++
	./gcov_report
	lcov -t "stest" -o s21_test.info -c -d . --ignore-errors mismatch
	genhtml -o report s21_test.info
	open ./report/index.html
	
clean:
	@-rm -rf *.o *.gcno *.gcda *.gcov *.info coverage_report *.a test gcov_report -r test.dSYM -r report

style:
	@clang-format -style=Google -n *.cpp *.h

leaks: $(LIB).a
	$(CC) $(FLAGS) -fsanitize=address -g -c $(FILES)
	ar rcs $(LIB).a *.o
	$(CC) -fsanitize=address -g $(FLAGS) $(LIB).a test.cpp -o test $(TEST_LIB)
	leaks -atExit -- ./test
