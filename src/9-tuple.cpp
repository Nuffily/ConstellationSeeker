#pragma warning(disable : 4146 4244 4828)

#include <iostream>
#include <vector>
#include <thread>
#include <mutex>
#include <gmpxx.h>
#include <chrono>

using namespace std;

mutex cout_mutex;

auto commence = chrono::high_resolution_clock::now();

// Возвращает вектор простых чисел в промежутке от start до limit
const vector<short> sieve(int limit, int start) {
    vector<bool> is_prime(limit + 1, true);
    vector<short> primes;
    for (int p = 2; p * p <= limit; ++p) {
        if (is_prime[p]) {
            for (int i = p * p; i <= limit; i += p)
                is_prime[i] = false;
        }
    }
    for (short p = start; p <= limit; ++p)
        if (is_prime[p])
            primes.push_back(p);
    return primes;
}

// Проверяет, делится ли n на одно из чисел в векторе
bool is_candidate(const mpz_class& n, const vector<short>& sie) {
    for (short p : sie) {
        if (mpz_divisible_ui_p(n.get_mpz_t(), p))
            return false;
    }
    return true;
}

// Проверяет, делится ли n на одно из чисел в векторе
bool is_candidate(const unsigned long& n, const vector<short>& sie) {
    for (short p : sie) {
        if (!(n % p)) {
            return false;
        }
    }
    return true;
}

/*
 Ищет девятки простых чисел, начиная со start, выводит всю девятку и время, прошедшее с запуска программы
    hard_sieve - вектор простостых чисел, с помощью которых будет проверяться простота чисел
    
Поиск проходит только по особым числам с особым шагом
    step - диаметр одного цикла поиска
    offs - промежутки между числами, которые нужно проверять (с которых может начинаться девятка)
    offs_last - offs.size - 1
    start_off - длина от начала цикла до первого нужного числа
    const_offset - показыват шаблон девятки, то есть какая длина между числами в девятке
 */
void search_quadruplets(const mpz_class& start, int thread_id, const vector<short>& hard_sieve,
    const short start_off, const short step, const vector<short>& offs, const int offs_last, const std::array<int,8>& const_offset) {
    
    mpz_class current[9];
    current[0] = start + (step - start % step) + start_off - step;

    unsigned int prime_mul = 1;

    const vector<short> light_sieve = sieve(40, 14);

    for (auto& prime : light_sieve)
        prime_mul *= prime;

    unsigned long checker = (prime_mul < current[0]) ? mpz_class(current[0] % prime_mul).get_si() : current[0].get_si();

    short iteration = 0;


    while (true) {

        if (is_candidate(checker, light_sieve)) {

            for (int i = 1; i <= 8; i++) {
                if (!is_candidate(checker + const_offset[i - 1], light_sieve)) {
                    goto exit_l;
                }
            }

            if (!is_candidate(current[0], hard_sieve)) {
                goto exit_l;
            }

            for (int i = 1; i <= 8; i++) {
                current[i] = current[0] + const_offset[i - 1];

                if (!is_candidate(current[i], hard_sieve)) {
                    goto exit_l;
                }
            }

            for (int i = 0; i <= 8; i++) {

                if (!mpz_probab_prime_p(current[i].get_mpz_t(), 1)) {
                    goto exit_l;
                }
            }

            for (int i = 0; i <= 8; i++) {

                if (!mpz_probab_prime_p(current[i].get_mpz_t(), 24))
                    goto exit_l;
            }

            {
            lock_guard<mutex> lock(cout_mutex);
                
            if (current[0] < start)
                goto exit_l;
                
            cout << "Thread " << thread_id << " found:\n";
            for (int i = 0; i < 9; i++)
                cout << current[i] << "\n";

            std::chrono::duration<double> dur = std::chrono::high_resolution_clock::now() - commence;
            cout << "\n" << dur.count() << "\n\n";
            }
            
        }
        exit_l:

        current[0] += offs[iteration];
        checker += offs[iteration];

        if (iteration < offs_last)
            iteration++;
        else
            iteration = 0;

        if (checker > prime_mul)
            checker -= prime_mul;
    }
}

int main() {
    mpz_class base = 10;
    unsigned long exp = 15;

    mpz_class z;
    mpz_pow_ui(z.get_mpz_t(), base.get_mpz_t(), exp);


    const vector<short> super_sie = sieve(400, 40);

    
    const vector<short> offs_1 = { 2730, 1890, 2730, 840, 3360, 1260, 1050, 420, 2310, 1890, 1260, 1470, 840, 4620, 3360};
    const std::array<int,8> const_offfs_1 = { 2, 6, 8, 12, 18, 20, 26, 30 };
    const short start_1 = 1271;
    const short step_1 = 30030;

    const vector<short> offs_2 = { 990, 420, 1320, 630, 360, 900, 990, 60, 1680, 990, 630, 1320, 1260, 2040, 270, 420, 2310, 570, 1320, 1410, 630, 690, 570, 60, 990, 1680, 1950, 990, 1260, 1320 };
    const std::array<int,8> const_offfs_2 = { 4, 6, 10, 16, 18, 24, 28, 30 };
    const short start_2 = 1273;
    const short step_2 = 30030;

    const vector<short> offs_3 = { 1260, 990, 1950, 1680, 990, 60, 570, 690, 630, 1410, 1320, 570, 2310, 420, 270, 2040, 1260, 1320, 630, 990, 1680, 60, 990, 900, 360, 630, 1320, 420, 990, 1320 };
    const std::array<int,8> const_offfs_3 = { 2, 6, 12, 14, 20, 24, 26, 30 };
    const short start_3 = 17;
    const short step_3 = 30030;

    const vector<short> offs_4 = { 4620, 840, 1470, 1260, 1890, 2310, 420, 1050, 1260, 3360, 840, 2730, 1890, 2730, 3360 };
    const std::array<int,8> const_offfs_4 = { 4, 10, 12, 18, 22, 24, 28, 30 };
    const short start_4 = 2059;
    const short step_4 = 30030;

    thread t1([&]() { search_quadruplets(z, 1, super_sie, start_1, step_1, offs_1, 14, const_offfs_1); });
    thread t2([&]() { search_quadruplets(z, 2, super_sie, start_2, step_2, offs_2, 29, const_offfs_2); });
    thread t3([&]() { search_quadruplets(z, 3, super_sie, start_3, step_3, offs_3, 29, const_offfs_3); });
    thread t4([&]() { search_quadruplets(z, 4, super_sie, start_4, step_4, offs_4, 14, const_offfs_4); });
    
    t1.join();
    t2.join();
    t3.join();
    t4.join();

    return 0;
}
