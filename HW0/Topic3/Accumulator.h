#pragma once

class Accumulator {
	float sum;
public:
	Accumulator() { sum = 0; }
	~Accumulator() {}

	float operator()(float to_add)
	{
		this->sum += to_add;
		return sum;
	}
};


