#pragma once

#include <random>;
#include <ctime>;

#include "boost/lexical_cast.hpp";


struct randomSequence
{
	uint64_t size;

	uint8_t* array;
	
	randomSequence(const uint64_t& size);
	uint8_t* getBitSequence();

	~randomSequence();
};
