#pragma once

#include <iostream>;
#include <vector>;
#include <algorithm>;
#include <cmath>;
#include <array>;
#include <string>;
#include <map>;
#include <numeric>;

#include <boost/math/distributions/normal.hpp>;

#include "./utils.hpp";
#include "./fast-fourier-transform.hpp";


static boost::math::normal dist(0.0, 1.0);


std::map<std::string, std::vector<double>> getChiSquareParameters(const uint64_t& l, const std::vector<double>& significanceLevels);

std::vector<double> equiprobabilityTest(const randomSequence* sequenceItem, const std::vector<double>& significanceLevels);
std::vector<double> independenceTest(const randomSequence* sequenceItem, const std::vector<double>& significanceLevels);
std::vector<double> uniformityTest(const randomSequence* sequenceItem, const uint32_t& segmentsCount, std::vector<double>& significanceLevels);

std::map<std::string, std::vector<double>> getSpectralTestParameters(const uint64_t& bitSize, const std::vector<double>& significanceLevels);
std::vector<double> spectralTest(randomSequence* sequenceItem, const std::vector<double>& significanceLevels);

uint8_t getL(const uint64_t& bitSize);
std::map<std::string, std::vector<double>> getUniversalTestParameters(const uint64_t& bitSize, const std::vector<double>& significanceLevels);
std::vector<double> universalTest(randomSequence* sequenceItem, const std::vector<double>& significanceLevels);

std::vector<double> processNormalDistributionCriteria(const std::vector<double>& statisticValues, const std::vector<double>& means, const std::vector<double>& dispersions, const std::vector<double>& significanceLevels);

void showTestResult(const std::vector<double>& statisticValues, const std::vector<double>& limitValues, const std::vector<double>& significanceLevels);
