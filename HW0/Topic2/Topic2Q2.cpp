#include <vector>
#include <iostream>
#include <algorithm>

int abs_value(int val)
{
	return (val < 0) ? -val : val;
}

int main()
{
	auto abs_compare = [](int val1, int val2)
	{
		return (abs_value(val1) < abs_value(val2)) ? true : false;
	};
	int arr[10] = { 1, 4, 2, -321, -3, 5, 8, 90, -89, 0 };
	std::cout << "The original array is - ";
	for (auto i : arr)
		std::cout << i << ' ';

	// Sort
	std::sort(arr, arr + 10, abs_compare);

	std::cout << std::endl << "The sorted array is - ";
	for (auto i : arr)
		std::cout << i << ' ';

	return 0;
}
