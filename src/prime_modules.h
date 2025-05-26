#pragma warning(disable : 4146 4244 4828)
#pragma once

#include <vector>
#include <gmpxx.h>

using namespace std;

const vector<int> sieve(int limit, int from);

// Проверка на делимость на числа из sieve
bool is_candidate(const mpz_t n, const vector<int>& sieve);

// Проверка на делимость на числа из sieve
bool is_candidate(const unsigned int& n, const vector<int>& sieve);

bool is_candidate(const mpz_class& n, const vector<int>& sie);

// Strong Prime Probability Test 
bool SPP2(mpz_t n);
