#pragma once
#pragma warning(disable : 4146 4244 4828)

#include <vector>
#include <gmpxx.h>


using namespace std;

void search_quadruplets(const mpz_class& start, int thread_id, const vector<int>& hard_sieve);

int seekForQuadruplets(mpz_class z);
