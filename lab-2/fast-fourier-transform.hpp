#pragma once

#include <cmath>;
#include <memory>;

#include "utils.hpp";


std::shared_ptr<double[]> executeFastFourierTransform(randomSequence* sequence, const uint64_t& bitSize);
