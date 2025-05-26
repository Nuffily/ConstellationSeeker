#pragma warning(disable : 4146 4244 4828)
#pragma once

#include <iostream>
#include <vector>
#include <thread>
#include <mutex>
#include <gmpxx.h>
#include <chrono>
#include <array>
#include <set>
#include <atomic>
#include <csignal>
#include <windows.h>
#include <cstdlib>
#include <cmath>
#include "prime_modules.h"

using namespace std;

mutex cout_mutex_4;

chrono::high_resolution_clock::time_point commence_4;

/*
    Поля для остановки программы
    При нажатии CTRL + C после некоторой задержки выводится последнее пройденное значение, на котором завершилась программа
*/
bool exit_flag_4 = false;
mpz_class min_checked_4("0");

void signal_handler_4(int signal) {
    exit_flag_4 = true;
}

void little_quadruplets(const mpz_class& start) {
    mpz_class current = start;
    const int first_quadruplets[12] = { 5, 11, 101, 191, 821, 1481, 1871, 2081, 3251, 3461, 5651, 9431 };

    for (int first_quaruplet : first_quadruplets) {
        if (current < first_quaruplet) {
            lock_guard<mutex> lock(cout_mutex_4);

            std::chrono::duration<double> dur = std::chrono::high_resolution_clock::now() - commence_4;

            cout << dur.count() << " : Thread " << 1 << " found:\n"
                << first_quaruplet << "\n" << first_quaruplet + 2 << "\n"
                << first_quaruplet + 6 << "\n" << first_quaruplet + 8 << "\n\n";

        }
    }
}

/*
 Ищем четверки простых чисел
 В каждом промежутке диаметром 2310 смотрим только числа с определенным шагом,
 остальные не допустимы, так как хоть один элемент четверки делился бы на 2, 3, 5, 7 или 11

 cheker - число, равное остатку current (текущего рассматриваемого числа) на prime_mul - произведение некоторых простых чисел (сейчас - от 13 до 36)

 cheker - переменная типа int, поэтому проверять ее делимость на малые простые быстрее, чем самого current, а их делимость на заданные простые совпадает

 Если все checker { +0, +2, +6, +8} числа прошли делимость на малые простые, начинаем поочередно делить числа предпологаемой четверки
    на простые числа из hard_sieve

 Затем применяем для них СПП2 тест

 Если все числа вновь прошли проверку, поочередно проводим для них тест Миллера-Рабина по одному разу

 Если и это испытание они прошли, скорее всего, они простые, но для уверенности для каждого запускае еще 24 теста Миллера-Рабина

 */

void search_quadruplets(const mpz_class& start, int thread_id, const vector<int>& hard_sieve) {

    const unsigned short next[21] = { 90, 30, 210, 90, 90, 210, 30, 90, 90, 120, 120, 90, 90, 30, 210, 90, 90, 210, 30, 90, 11760 };
    int iteration = 0;
    mpz_class temp;

    if (start < 10000) {
        if (thread_id == 1) // При вводе слишком маленького числа выводим младшие четверки через little_quadruplets, делать это должен лишь один поток. А поиск начинаем с 10000
            little_quadruplets(start);

        temp = 10000 - (10000 % 2310) + 101 + 2310 * (thread_id - 7);
    }
    else {
        temp = start - (start % 2310) + 101 + 2310 * (thread_id - 7);
    }

    mpz_t current, current2, current8, current6;


    mpz_init(current);
    mpz_init(current2);
    mpz_init(current6);
    mpz_init(current8);

    mpz_set(current, temp.get_mpz_t());

    unsigned int prime_mul = 1;
    unsigned int billion = 0;

    const vector<int> light_sieve = sieve(36, 13);

    for (auto& prime : light_sieve)
        prime_mul *= prime;

    const unsigned int thershold = 1041872676;

    unsigned int checker = mpz_class(temp % prime_mul).get_si();


    while (true) {

        mpz_add_ui(current, current, next[iteration]);
        checker += next[iteration];

        if (iteration > 19) {

            iteration = 0;

            if (checker > thershold) {
                checker -= thershold;

                if (thread_id == 1) {
                    billion++;            // При каждом пройденном миллиарде чисел выводим сообщение об этом в консоль
                    std::chrono::duration<double> dur = std::chrono::high_resolution_clock::now() - commence_4;
                    cout << dur.count() << " : " << billion << " Billion is passed\n" << endl;
                }
            }
        }
        else {
            iteration++;
            if (exit_flag_4) {
                lock_guard<mutex> lock(cout_mutex_4);

                if (min_checked_4 == 0 || min_checked_4 > mpz_class(current))
                    min_checked_4 = mpz_class(current);

                return;
            }
        }

        // Делим на малые простые
        if (!is_candidate(checker, light_sieve)) continue;

        if (!is_candidate(checker + 8, light_sieve))  continue;

        if (!is_candidate(checker + 6, light_sieve))  continue;

        if (!is_candidate(checker + 2, light_sieve))  continue;


        // Делим на не очень малые простые
        if (!is_candidate(current, hard_sieve)) continue;

        mpz_add_ui(current8, current, 8);
        if (!is_candidate(current8, hard_sieve)) continue;

        mpz_add_ui(current6, current, 6);
        if (!is_candidate(current6, hard_sieve)) continue;

        mpz_add_ui(current2, current, 2);
        if (!is_candidate(current2, hard_sieve)) continue;

        // Тест СПП2
        if (!(SPP2(current) && SPP2(current8) && SPP2(current6) && SPP2(current2))) continue;

        // Тест Миллера-Рабина
        if (!mpz_probab_prime_p(current, 1)) continue;

        if (!mpz_probab_prime_p(current8, 1)) continue;

        if (!mpz_probab_prime_p(current6, 1)) continue;

        if (!mpz_probab_prime_p(current2, 1)) continue;

        if (mpz_probab_prime_p(current6, 24) &&
            mpz_probab_prime_p(current8, 24) &&
            mpz_probab_prime_p(current2, 24) &&
            mpz_probab_prime_p(current, 24)) {

            if (mpz_cmp(start.get_mpz_t(), current)) {
                lock_guard<mutex> lock(cout_mutex_4);

                std::chrono::duration<double> dur = std::chrono::high_resolution_clock::now() - commence_4;

                cout << dur.count() << " : Thread " << thread_id << " found:\n"
                    << current << "\n" << current2 << "\n"
                    << current6 << "\n" << current8 << "\n\n";
            }
        }

    }
}



int seekForQuadruplets(mpz_class z) {

    signal(SIGINT, signal_handler_4);                         // Для перехвата CTRL + C

    //mpz_class base = 10;
    //unsigned long exp = 1000;

    //mpz_class z;
    //mpz_pow_ui(z.get_mpz_t(), base.get_mpz_t(), exp);       // Текущее число = base ^ exp

    // z += mpz_class("140358297211");                        // Добавка (удобно для продолжения с последней точки)

    int digits = mpz_sizeinbase(z.get_mpz_t(), 10);

    if (digits < 10)
        digits = 10;

    commence_4 = std::chrono::high_resolution_clock::now();

    const vector<int> hard_sieve = sieve(digits * sqrt(digits) / 1.75, 37);  // Вектор с простыми числами для проверки

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

    cout << "Program has stopped on number\n" << min_checked_4 << endl;            // Информация о последнем пройденном числе

    return 0;
}
