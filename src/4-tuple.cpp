
#pragma warning(disable : 4146 4244 4828)

#include <iostream>
#include <vector>
#include <thread>
#include <mutex>
#include <gmpxx.h>
#include <chrono>

using namespace std;

mutex cout_mutex;

auto commence = std::chrono::high_resolution_clock::now();

// Возвращает вектор простых чисел в промежутке от from до limit
const vector<int> sieve(int limit,int from) {
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
bool is_candidate(const mpz_class& n, const vector<int>& sieve) {
    for (int p : sieve) {
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

/*
 Ищем четверки простых чисел
 В каждом промежутке диаметром 2310 смотрим только числа с определенным шагом,
 остальные не допустимы, так как хоть один элемент четверки делился бы на 2, 3, 5, 7 или 11
 
 cheker - число, равное остатку current (текущего рассматриваемого числа) на prime_mul - произведение некоторых простых чисел (сейчас - от 13 до 36)
 
 cheker - переменная типа int, поэтому проверять ее делимость на малые простые быстрее, чем самого current, а их делимость на заданные простые совпадает
 
 Если все checker { +0, +2, +6, +8} числа прошли делимость на малые простые, начинаем поочередно делить числа предпологаемой четверки
    на простые числа из hard_sieve
 
 Если все числа вновь прошли проверку, поочередно проводим для них тест Миллера-Рабина по одному разу
 
 Если и это испытание они прошли, скорее всего, они простые, но для уверенности для каждого запускае еще 24 теста Миллера-Рабина

 */

void search_quadruplets(const mpz_class& start, int thread_id, const vector<int>& hard_sieve) {

    const short next[21] = { 90, 30, 210, 90, 90, 210, 30, 90, 90, 120, 120, 90, 90, 30, 210, 90, 90, 210, 30, 90, 11760 };
    int iteration = 0;

    mpz_class current = start - (start % 2310) + 101 + 2310 * (thread_id - 1);
    mpz_class current2, current8, current6;

    unsigned int prime_mul = 1;

    const vector<int> light_sieve = sieve(36, 13);

    for (auto& prime : light_sieve)
        prime_mul *= prime;

    unsigned int checker = mpz_class(current % prime_mul).get_si();

    while (true) {

        if (is_candidate(checker, light_sieve)) {

            if (is_candidate(checker + 8, light_sieve)) {
                
                if (is_candidate(checker + 6, light_sieve)) {
                    
                    if (is_candidate(checker + 2, light_sieve)) {

                        if (is_candidate(current, hard_sieve)) {
                            current8 = current + 8;
                            if (is_candidate(current8, hard_sieve)) {
                                current6 = current + 6;
                                if (is_candidate(current6, hard_sieve)) {
                                    current2 = current + 2;
                                    if (is_candidate(current2, hard_sieve)) {

                                        if (mpz_probab_prime_p(current.get_mpz_t(), 1)) {

                                            if (mpz_probab_prime_p(current8.get_mpz_t(), 1)) {

                                                if (mpz_probab_prime_p(current6.get_mpz_t(), 1)) {

                                                    if (mpz_probab_prime_p(current2.get_mpz_t(), 1)) {

                                                        if (mpz_probab_prime_p(current6.get_mpz_t(), 24) &&
                                                            mpz_probab_prime_p(current8.get_mpz_t(), 24) &&
                                                            mpz_probab_prime_p(current2.get_mpz_t(), 24) &&
                                                            mpz_probab_prime_p(current.get_mpz_t(), 24)) {

                                                            if (current > start) {
                                                                lock_guard<mutex> lock(cout_mutex);
                                                                
                                                                cout << "Thread " << thread_id << " found:\n"
                                                                << current << "\n" << current + 2 << "\n"
                                                                << current + 6 << "\n" << current + 8 << "\n\n";
                                                                
                                                                std::chrono::duration<double> dur = std::chrono::high_resolution_clock::now() - commence;
                                                                
                                                                cout << dur.count() << "\n\n";
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        current += next[iteration];
        checker += next[iteration];
        
        if (iteration > 19)
            iteration = 0;
        else
            iteration++;

        if (checker > prime_mul)
            checker -= prime_mul;
    }
}

int main() {
    
    mpz_class base = 13;
    unsigned long exp = 500;

    mpz_class z;
    mpz_pow_ui(z.get_mpz_t(), base.get_mpz_t(), exp);       // Текущее число = base ^ exp

    const vector<int> hard_sieve = sieve(1000, 36);
    
    thread t1([&]() { search_quadruplets(z, 1, hard_sieve); });
    thread t2([&]() { search_quadruplets(z, 2, hard_sieve); });
    thread t3([&]() { search_quadruplets(z, 3, hard_sieve); });
    thread t4([&]() { search_quadruplets(z, 4, hard_sieve); });
    thread t5([&]() { search_quadruplets(z, 5, hard_sieve); });
    thread t6([&]() { search_quadruplets(z, 6, hard_sieve); });

    t1.join();
    t2.join();
    t3.join();
    t4.join();
    t5.join();
    t6.join();

    return 0;
}
