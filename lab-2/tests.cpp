#include "tests.hpp";


void showTestResult(const std::vector<double>& statisticValues, const std::vector<double>& limitValues, const std::vector<double>& significanceLevels)
{
	for (uint8_t i = 0; i < statisticValues.size(); ++i)
	{
		auto significanceLevel = significanceLevels[i];
		auto statisticValue = statisticValues[i];
		auto limitValue = limitValues[i];

		auto predicateValue = statisticValue <= limitValue;
		
		printf("Significance level = %f, statistic value = %f, limit value = %f, test passed = %s\n", significanceLevel, statisticValue, limitValue, predicateValue ? "\033[32mtrue\033[0m" : "\033[31mfalse\033[0m");
	}
}

std::map<std::string, std::vector<double>> getChiSquareParameters(const uint64_t& l, const std::vector<double>& significanceLevels)
{
	uint8_t testCasesCount = significanceLevels.size();

	auto limitValues = std::vector<double>(testCasesCount);

	for (uint8_t i = 0; i < testCasesCount; ++i)
	{
		limitValues[i] = sqrt(2 * l) * quantile(dist, 1 - significanceLevels[i]) + l;
	}

	return std::map<std::string, std::vector<double>> { { "limits", limitValues } };
}

std::vector<double> equiprobabilityTest(const randomSequence* sequenceItem, const std::vector<double>& significanceLevels)
{
	auto byteSize = sequenceItem->size >> 3;
	auto byteSequence = sequenceItem->array;

	auto bytesFrequency = std::array<int64_t, UCHAR_MAX + 1>();

	int64_t expectation = byteSize / 256;

	for (uint64_t i = 0; i < byteSize; ++i)
	{
		++bytesFrequency[byteSequence[i]];
	}

	double chi_2 = 0;

	for (uint64_t i = 0; i < bytesFrequency.size(); ++i)
	{
		chi_2 += pow(bytesFrequency[i] - expectation, 2) / expectation;
	}

	return std::vector<double>(significanceLevels.size(), chi_2);
}

std::vector<double> independenceTest(const randomSequence* sequenceItem, const std::vector<double>& significanceLevels)
{
	auto bitSize = sequenceItem->size;
	auto byteSequence = sequenceItem->array;

	auto bytesPairs = std::array<std::array<uint64_t, UCHAR_MAX + 1>, UCHAR_MAX + 1>();
	auto byteFirstPos = std::array<uint64_t, UCHAR_MAX + 1>();
	auto byteSecondPos = std::array<uint64_t, UCHAR_MAX + 1>();

	for (uint64_t i = 0; i < bitSize >> 3; i += 2)
	{
		++bytesPairs[byteSequence[i]][byteSequence[i + 1]];
		++byteFirstPos[byteSequence[i]];
		++byteSecondPos[byteSequence[i + 1]];
	}

	double sum = 0;

	for (uint64_t i = 0; i <= UCHAR_MAX; ++i)
	{
		for (uint64_t j = 0; j <= UCHAR_MAX; ++j)
		{
			if (byteFirstPos[i] * byteSecondPos[j] != 0)
			{
				sum += pow(bytesPairs[i][j], 2) / (byteFirstPos[i] * byteSecondPos[j]);
			}
		}
	}

	double chi_2 = (bitSize / 16) * (sum - 1);

	return std::vector<double>(significanceLevels.size(), chi_2);
}

std::vector<double> uniformityTest(const randomSequence* sequenceItem, const uint32_t& segmentsCount, std::vector<double>& significanceLevels)
{
	auto byteSize = sequenceItem->size >> 3;
	auto byteSequence = sequenceItem->array;

	auto bytesFrequency = std::vector<std::array<uint64_t, UCHAR_MAX + 1>>(segmentsCount);
	for (uint32_t i = 0; i < segmentsCount; ++i)
	{
		bytesFrequency[i] = std::array<uint64_t, UCHAR_MAX + 1>();
	}

	uint32_t segmentSize = byteSize / segmentsCount;

	for (uint64_t i = 0; i < segmentsCount; ++i)
	{
		for (uint32_t j = 0; j < segmentSize; ++j)
		{
			++bytesFrequency[i][byteSequence[i * segmentSize + j]];
		}
	}

	auto bytesCommonFrequency = std::array<uint64_t, UCHAR_MAX + 1>();

	for (uint64_t i = 0; i < segmentsCount; ++i)
	{
		for (uint32_t j = 0; j <= UCHAR_MAX; ++j)
		{
			bytesCommonFrequency[j] += bytesFrequency[i][j];
		}
	}

	double sum = 0;

	for (uint64_t i = 0; i < segmentsCount; ++i)
	{
		for (uint32_t j = 0; j <= UCHAR_MAX; ++j)
		{
			if (bytesCommonFrequency[j] != 0) sum += pow(bytesFrequency[i][j], 2) / static_cast<double>(bytesCommonFrequency[j] * segmentSize);
		}
	}

	auto chi_2 = byteSize * (sum - 1);

	return std::vector<double>(significanceLevels.size(), chi_2);
}

std::map<std::string, std::vector<double>> getSpectralTestParameters(const uint64_t& bitSize, const std::vector<double>& significanceLevels)
{
	uint8_t testCasesCount = significanceLevels.size();

	auto means = std::vector<double>(testCasesCount);
	auto dispersions = std::vector<double>(testCasesCount);
	auto limitValues = std::vector<double>(testCasesCount);

	for (uint8_t i = 0; i < testCasesCount; ++i)
	{
		auto significanceLevel = significanceLevels[i];

		means[i] = (1 - significanceLevel) * bitSize / 2.0;
		dispersions[i] = sqrt(bitSize * 0.25 * (1 - significanceLevel) * significanceLevel);
		limitValues[i] = quantile(dist, 1 - significanceLevels[i]);
	}

	return std::map<std::string, std::vector<double>> { {"means", means}, { "dispersions", dispersions }, { "limits", limitValues } };
}

std::vector<double> spectralTest(randomSequence* sequenceItem, const std::vector<double>& significanceLevels)
{
	auto bitSize = sequenceItem->size;

	auto X = executeFastFourierTransform(sequenceItem, bitSize);

	auto m = std::vector<double>(bitSize / 2 + 1);
	m[0] = std::abs(X[0]);

	for (uint64_t i = 1; i <= bitSize / 2; ++i)
	{
		m[i] = sqrt(pow(X[2 * i - 1], 2) + pow(X[2 * i], 2));
	}

	uint8_t testCasesCount = significanceLevels.size();

	auto statisticValues = std::vector<double>(testCasesCount);

	for (uint8_t i = 0; i < testCasesCount; ++i)
	{
		auto upperBound = sqrt(log(1 / significanceLevels[i]) * bitSize);

		statisticValues[i] = std::accumulate(m.begin(), m.end(), 0, [&](uint64_t accumulator, double currentItem) {
			return accumulator + (currentItem < upperBound ? 1 : 0);
		});
	}

	return statisticValues;
}

uint8_t getL (const uint64_t& bitSize)
{
	if (bitSize >= 1059061760)	return 16;
	if (bitSize >= 496435200)	return 15;
	if (bitSize >= 231669760)	return 14;
	if (bitSize >= 107560960)	return 13;
	if (bitSize >= 49643520)	return 12;
	if (bitSize >= 22753280)	return 11;
	if (bitSize >= 10342400)	return 10;
	if (bitSize >= 4654080)		return 9;
	if (bitSize >= 2068480)		return 8;
	if (bitSize >= 904960)		return 7;
	if (bitSize >= 387840)		return 6;
	
	return 5;
}

std::map<std::string, std::vector<double>> getUniversalTestParameters(const uint64_t& bitSize, const std::vector<double>& significanceLevels)
{
	double	expected_value[17] = {
		0, 0, 0, 0, 0, 0, 5.2177052, 6.1962507, 7.1836656, 8.1764248, 9.1723243, 10.170032, 11.168765, 12.168070, 13.167693, 14.167488, 15.167379
	};
	double variance[17] = {
		0, 0, 0, 0, 0, 0, 2.954, 3.125, 3.238, 3.311, 3.356, 3.384, 3.401, 3.410, 3.416, 3.419, 3.421
	};

	auto L = getL(bitSize);
	auto p = pow(2, L);

	auto Q = 10 * p;
	auto K = floor(bitSize / L) - Q;

	auto c = 0.7 - 0.8 / L + (4 + 32.0 / L) * pow(K, -3.0 / L) / 15;

	uint8_t testCasesCount = significanceLevels.size();

	auto means = std::vector<double>(testCasesCount, expected_value[L]);
	auto dispersions = std::vector<double>(testCasesCount, c * sqrt(variance[L] / K));
	auto limitValues = std::vector<double>(testCasesCount);
	
	for (uint8_t i = 0; i < testCasesCount; ++i)
	{
		limitValues[i] = quantile(dist, 1 - significanceLevels[i]);
	}

	return std::map<std::string, std::vector<double>> { {"means", means}, { "dispersions", dispersions }, { "limits", limitValues } };
}

std::vector<double> universalTest(randomSequence* sequenceItem, const std::vector<double>& significanceLevels)
{
	auto bitSize = sequenceItem->size;
	auto bitSequence = std::unique_ptr<uint8_t[]>(sequenceItem->getBitSequence());

	auto L = getL(bitSize);
	auto p = pow(2, L);

	auto Q = 10 * p;
	auto K = floor(bitSize / L) - Q;

	auto T = std::vector<int64_t>(p, 0);

	auto sum = 0.0;

	for (uint64_t i = 1; i < Q; ++i)
	{
		double decRep = 0;

		for (uint64_t j = 0; j < L; ++j)
		{
			decRep += bitSequence[(i - 1) * L + j] * pow(2, L - 1 - j);
		}

		T[decRep] = i;
	}

	for (uint64_t i = Q + 1; i <= Q + K; ++i) {
		double decRep = 0;

		for (uint64_t j = 0; j < L; ++j)
		{
			decRep += bitSequence[(i - 1) * L + j] * pow(2, L - 1 - j);
		}

		sum += log(i - T[decRep]) / log(2);
		
		T[decRep] = i;
	}

	return std::vector<double>(significanceLevels.size(), sum / K);
}

std::vector<double> processNormalDistributionCriteria(const std::vector<double>& statisticValues, const std::vector<double>& means, const std::vector<double>& dispersions, const std::vector<double>& significanceLevels)
{
	uint8_t testCasesCount = significanceLevels.size();

	auto normalizedStatisticValue = std::vector<double>(testCasesCount);

	for (uint8_t i = 0; i < testCasesCount; ++i)
	{
		normalizedStatisticValue[i] = std::abs(statisticValues[i] - means[i]) / dispersions[i];
	}

	return normalizedStatisticValue;
}
