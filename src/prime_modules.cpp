#pragma warning(disable : 4146 4244 4828)

#include <iostream>
#include <vector>
#include <gmpxx.h>

using namespace std;
// Возвращает вектор простых чисел в промежутке от from до limit
const vector<int> sieve(int limit, int from) {
    vector<bool> is_prime(limit + 1, true);
    vector<int> primes;
    for (short p = 2; p * p <= limit; ++p) {
        if (is_prime[p]) {
            for (int i = p * p; i <= limit; i += p)
                is_prime[i] = false;
        }
    }
    for (int p = from; p <= limit; ++p)
        if (is_prime[p])
            primes.push_back(p);

    return primes;
}

// Проверка на делимость на числа из sieve
bool is_candidate(const mpz_t n, const vector<int>& sieve) {
    for (unsigned int p : sieve) {
        if (mpz_divisible_ui_p(n, p))
            return false;
    }
    return true;
}

// Проверка на делимость на числа из sieve
bool is_candidate(const mpz_class& n, const vector<int>& sie) {
    for (short p : sie) {
        if (mpz_divisible_ui_p(n.get_mpz_t(), p))
            return false;
    }
    return true;
}

// Проверка на делимость на числа из sieve
bool is_candidate(const unsigned int& n, const vector<int>& sieve) {
    for (int p : sieve) {
        if (!(n % p)) {
            return false;
        }
    }
    return true;
}

// Strong Prime Probability Test 
bool SPP2(mpz_t n) {
    mpz_t t;
    mpz_t d;
    mpz_t ff;
    mpz_t n_minus_1;
    mpz_t base_2;

    mpz_inits(t, d, ff, n_minus_1, base_2, NULL);

    mpz_sub_ui(n_minus_1, n, 1);

    mpz_set(t, n_minus_1);

    unsigned long s = mpz_scan1(t, 0);

    mpz_fdiv_q_2exp(d, t, s);

    mpz_set_ui(base_2, 2);
    mpz_powm(ff, base_2, d, n);

    if (mpz_cmp_ui(ff, 1) == 0 || mpz_cmp(ff, n_minus_1) == 0) {
        mpz_clears(t, d, ff, n_minus_1, base_2, NULL);
        return true;
    }

    for (unsigned long i = 1; i < s; ++i) {
        mpz_powm_ui(ff, ff, 2, n);
        if (mpz_cmp(ff, n_minus_1) == 0) {
            mpz_clears(t, d, ff, n_minus_1, base_2, NULL);
            return true;
        }
    }

    mpz_clears(t, d, ff, n_minus_1, base_2, NULL);
    return false;
}
