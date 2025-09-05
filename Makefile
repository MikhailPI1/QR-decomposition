1CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++11 -O2 -Iinclude

# Директории
SRC_DIR = src
APP_DIR = apps
TEST_DIR = tests
INCLUDE_DIR = include
BUILD_DIR = build


MAIN_TARGET = $(BUILD_DIR)/qr_program.exe


UNIT_TEST_TARGET = $(BUILD_DIR)/unit_tests.exe
INTEGRATION_TEST_TARGET = $(BUILD_DIR)/integration_tests.exe

HEADERS = $(INCLUDE_DIR)/qr_decomposition.h


$(shell if not exist $(BUILD_DIR) mkdir $(BUILD_DIR))


all: $(MAIN_TARGET)

test: $(UNIT_TEST_TARGET) $(INTEGRATION_TEST_TARGET)
	./$(UNIT_TEST_TARGET)
	./$(INTEGRATION_TEST_TARGET)

test-unit: $(UNIT_TEST_TARGET)
	./$(UNIT_TEST_TARGET)

test-integration: $(INTEGRATION_TEST_TARGET)
	./$(INTEGRATION_TEST_TARGET)


$(MAIN_TARGET): $(SRC_DIR)/qr_decomposition.cpp $(APP_DIR)/main.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $(SRC_DIR)/qr_decomposition.cpp $(APP_DIR)/main.cpp


$(UNIT_TEST_TARGET): $(SRC_DIR)/qr_decomposition.cpp $(TEST_DIR)/unit_tests.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $(SRC_DIR)/qr_decomposition.cpp $(TEST_DIR)/unit_tests.cpp


$(INTEGRATION_TEST_TARGET): $(SRC_DIR)/qr_decomposition.cpp $(TEST_DIR)/integration_tests.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $@ $(SRC_DIR)/qr_decomposition.cpp $(TEST_DIR)/integration_tests.cpp


clean:
	del /q /f $(BUILD_DIR)\* 2>nul
	if exist $(BUILD_DIR) rmdir /q /s $(BUILD_DIR)


rebuild: clean all


run: $(MAIN_TARGET)
	./$(MAIN_TARGET)

.PHONY: all clean rebuild run test test-unit test-integration