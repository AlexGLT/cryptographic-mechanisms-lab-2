#include <iostream>;
#include <fstream>;
#include <numeric>;

#include "utils.hpp";
#include "tests.hpp";
#include <openssl/rand.h>;


int main()
{
	std::vector<double> significanceLevels = { 0.01, 0.02, 0.05 };

	uint8_t levelsCount = significanceLevels.size();
	uint8_t iterationsCount = 10;
	uint64_t bitSize = static_cast<uint64_t>(1) << 30;

	std::vector<std::shared_ptr<randomSequence>> randomSequences(iterationsCount);

	for (uint8_t i = 0; i < iterationsCount; ++i)
	{
		randomSequences[i] = std::shared_ptr<randomSequence>(new randomSequence(bitSize));
		RAND_bytes(randomSequences[i]->array, bitSize >> 3);
	}

	printf("======OpenSSL Generator tests======\n");

	printf("1. Equiprobability test\n");

	auto equprobabilityTestParameters = getChiSquareParameters(UCHAR_MAX + 1, significanceLevels);

	auto equiprobabilityStatisticValues = std::vector<double>(levelsCount);

	for (uint8_t i = 0; i < iterationsCount; ++i)
	{
		auto statisticValues = equiprobabilityTest(randomSequences[i].get(), significanceLevels);

		for (uint8_t j = 0; j < levelsCount; ++j)
		{
			equiprobabilityStatisticValues[j] += statisticValues[j] / 10;
		}
	}

	showTestResult(equiprobabilityStatisticValues, equprobabilityTestParameters["limits"], significanceLevels);

	printf("\n");

	printf("2. Independence test\n");

	auto independenceTestParameters = getChiSquareParameters(pow(UCHAR_MAX + 1, 2), significanceLevels);

	auto independenceStatisticValues = std::vector<double>(levelsCount);

	for (uint8_t i = 0; i < iterationsCount; ++i)
	{
		auto statisticValues = independenceTest(randomSequences[i].get(), significanceLevels);

		for (uint8_t j = 0; j < levelsCount; ++j)
		{
			independenceStatisticValues[j] += statisticValues[j] / 10;
		}
	}

	showTestResult(independenceStatisticValues, independenceTestParameters["limits"], significanceLevels);

	printf("\n");

	printf("3. Uniformity test\n");

	auto segmentsCount = 16;
	auto uniformityTestParameters = getChiSquareParameters((UCHAR_MAX + 1) * (segmentsCount - 1), significanceLevels);

	auto uniformityStatisticValues = std::vector<double>(levelsCount);

	for (uint8_t i = 0; i < iterationsCount; ++i)
	{
		auto statisticValues = uniformityTest(randomSequences[i].get(), segmentsCount, significanceLevels);

		for (uint8_t j = 0; j < levelsCount; ++j)
		{
			uniformityStatisticValues[j] += statisticValues[j] / 10;
		}
	}

	showTestResult(uniformityStatisticValues, uniformityTestParameters["limits"], significanceLevels);

	printf("\n");

	printf("4. Spectral test\n");

	auto spectralTestParameters = getSpectralTestParameters(bitSize, significanceLevels);

	auto spectralStatisticValues = std::vector<double>(levelsCount);

	for (uint8_t i = 0; i < iterationsCount; ++i)
	{
		auto statisticValues = spectralTest(randomSequences[i].get(), significanceLevels);

		for (uint8_t j = 0; j < levelsCount; ++j)
		{
			spectralStatisticValues[j] += statisticValues[j] / 10;
		}
	}

	auto normalizedSpectralStatisticValue = processNormalDistributionCriteria(
		spectralStatisticValues, spectralTestParameters["means"], spectralTestParameters["dispersions"], significanceLevels
	);

	showTestResult(normalizedSpectralStatisticValue, spectralTestParameters["limits"], significanceLevels);

	printf("\n");

	printf("5. Universal test\n");

	auto universalTestParameters = getUniversalTestParameters(bitSize, significanceLevels);

	auto universalStatisticValues = std::vector<double>(levelsCount);

	for (uint8_t i = 0; i < iterationsCount; ++i)
	{
		auto statisticValues = universalTest(randomSequences[i].get(), significanceLevels);

		for (uint8_t j = 0; j < levelsCount; ++j)
		{
			universalStatisticValues[j] += statisticValues[j] / 10;
		}
	}

	auto normalizedUniversalStatisticValue = processNormalDistributionCriteria(
		universalStatisticValues, universalTestParameters["means"], universalTestParameters["dispersions"], significanceLevels
	);

	showTestResult(normalizedUniversalStatisticValue, universalTestParameters["limits"], significanceLevels);

	printf("\n");

	return 1;
}
