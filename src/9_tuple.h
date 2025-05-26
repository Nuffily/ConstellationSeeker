#pragma warning(disable : 4146 4244 4828)
#pragma once

#include <vector>
#include <gmpxx.h>

using namespace std;

void search_tuples_of_nine(const mpz_class& start, int thread_id, const vector<int>& hard_sieve,
    const int start_off, const int step, const vector<int>& offs, const int offs_last, const vector<int>& const_offset);

int seekForNineTuples(mpz_class z);
