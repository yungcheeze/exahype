
#include <iostream>

#include <chrono> // C++11

void chrono_example() {
	// this just uses C++11 std::chrono.
	std::cout << "Time since epoch, in periods: " << std::chrono::steady_clock::now().time_since_epoch().count() << "\n";
}
