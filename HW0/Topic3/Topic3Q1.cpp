#include "Accumulator.h"
#include <vector>
#include <cassert>
#include <iostream>

int main()
{
	// Q1.a
	std::vector<float> pdf_vector = { 0.1f,0.05f,0.4f,0.15f,0.2f,0.1f };

	// Q1.b
	float vec_sum = 0;
	for (float i : pdf_vector)
	{
		vec_sum += i;
	}
	assert(vec_sum == 1);

	// Q1.c
	Accumulator acc;
	std::vector<float> cdf_vector;
	for (float i : pdf_vector)
	{
		cdf_vector.push_back(acc(i));
	}
	
	// Prints
	std::cout << "The pdf vector is - ";
	for (float i : pdf_vector)
	{
		std::cout << i << " ";
	}

	std::cout << std::endl << "The cdf vector is - ";
	for (float i : cdf_vector)
	{
		std::cout << i << " ";
	}

	return 0;
}

//
//void perform_operation(int x, int y, int (*operation)(int, int)) {
//	std::cout << "Performing operation on " << x << " and " << y << ": ";
//	std::cout << operation(x, y) << std::endl;
//}
//





