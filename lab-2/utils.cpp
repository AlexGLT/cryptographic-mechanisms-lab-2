#include "utils.hpp";


randomSequence::randomSequence(const uint64_t& size)
{
    this->size = size;

    this->array = new uint8_t[this->size / 8];
}

uint8_t* randomSequence::getBitSequence()
{
    auto bitSequence = new uint8_t[this->size];

    for (uint64_t i = 0; i < this->size / 8; ++i)
    {
        for (uint8_t j = 0; j < 8; ++j)
        {
            bitSequence[8 * i + j] = (this->array[i] >> (7 - j)) & 1;
        }

    }

    return bitSequence;
}

randomSequence::~randomSequence()
{
    delete[] this->array;
}
