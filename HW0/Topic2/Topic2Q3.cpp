#include <vector>
#include <algorithm>
#include <iostream>

int main()
{
	int comparison_counter = 0;
	std::vector<int> vec{ 0,5,2,9,7,6,1,3,4,8 };
	auto sort_lambda = [&comparison_counter](int val1, int val2)
	{
		comparison_counter++;
		return (val1 < val2) ? true : false;
	};
	std::cout << "The original vector is - ";
	for (auto i : vec)
	{
		std::cout << i << ' ';
	}

	// sort
	std::sort(vec.begin(), vec.end(), sort_lambda);

	std::cout << std::endl << "The sorted vector is - ";
	for (auto i : vec)
	{
		std::cout << i << ' ';
	}

	std::cout << std::endl << "The number of comparisons is - " << comparison_counter << std::endl;
	return 0;
}