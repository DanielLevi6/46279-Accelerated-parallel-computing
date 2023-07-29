#include <memory>
#include <iostream>
#include "SmartVSRaw.h"

using namespace std;

void allocateRaw()
{
	int day = 20, month = 4, year = 2023;
	Date* rawPointer = new Date(day, month, year);
	cout << "The date today - " << rawPointer->getDay() << "/" << 
		rawPointer->getMonth() << "/" << rawPointer->getYear() << std::endl;

	delete rawPointer;
}

void allocateSmart()
{
	int day = 20, month = 4, year = 2023;
	unique_ptr<Date> smartDate(new Date(day, month, year));
	cout << "The date today - " << rawPointer->getDay() << "/" <<
		rawPointer->getMonth() << "/" << rawPointer->getYear() << std::endl;

}

int main()
{

	for (int i = 0; i < 1000000; i++)
	{
		allocateRaw();
	}
	for (int i = 0; i < 1000000; i++)
	{
		allocateSmart();
	}

	return 0;
}
