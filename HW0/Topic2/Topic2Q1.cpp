#include <vector>
#include <iostream>

int main()
{
	auto even_vector = [](std::vector<int>& nums_vector)
	{
		std::vector<int> new_vector;
		for (size_t i = 0; i < nums_vector.size(); i++)
		{
			if (nums_vector.at(i) % 2 == 0)
			{
				new_vector.insert(new_vector.end(), nums_vector.at(i));
			}
		}
		return new_vector;
	};
	std::vector<int> original_vector{ 13, 12, 41, 53, 64, 103, 132, 433213, 43232 };
	std::vector<int> new_vector = even_vector(original_vector);
	std::cout << "The original vector is - ";
	for (auto i : original_vector)
		std::cout << i << ' ';

	std::cout << std::endl << "The even vector is - ";
	for (auto i : new_vector)
		std::cout << i << ' ';

	return 0;
}