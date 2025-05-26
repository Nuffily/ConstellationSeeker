#pragma warning(disable : 4146 4244 4828)
#pragma once

#include <gmpxx.h>

#include "prime_modules.h"
#include "4_tuple.h"
#include "9_tuple.h"
#include <iostream>

using namespace std;

int main() {

	string input;
	mpz_class z;

	while (1) {

		cout << "What to seek? 4 for 4-tuples and 9 for 9-tuples: " << endl;

		cin >> input;

		if (input == "4") {

			cout << "Write a number in: " << endl;

			while (1) {

				cin >> input;

				try {
					z = mpz_class(input);
					break;
				}
				catch (exception) {
					cout << "Write a number in: " << endl;
				}
			}

			seekForQuadruplets(z);
			break;
		}

		else if (input == "9") {

			cout << "Write a number in: " << endl;

			while (1) {

				cin >> input;

				try {
					z = mpz_class(input);
					break;
				}
				catch (exception) {
					cout << "Write a number in: " << endl;
				}
			}

			seekForNineTuples(z);
			break;
		}	
	}

	return 0;
}
