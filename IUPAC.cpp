#include "IUPAC.hpp"

using namespace std;

const bool matrix[15][15] = {
	//         A      C      M      G      R      S      V      T      W      Y      H      K      D      B      N
	/* A */ {true , false, true , false, true , false, true , false, true , false, true , false, true , false, true }, // A
	/* C */ {false, true , true , false, false, true , true , false, false, true , true , false, false, true , true }, // C
	/* M */ {true , true , true , false, true , true , true , false, true , true , true , false, true , true , true }, // M
	/* G */ {false, false, false, true , true , true , true , false, false, false, false, true , true , true , true }, // G
	/* R */ {true , false, true , true , true , true , true , false, true , false, true , true , true , true , true }, // R
	/* S */ {false, true , true , true , true , true , true , false, false, true , true , true , true , true , true }, // S
	/* V */ {true , true , true , true , true , true , true , false, true , true , true , true , true , true , true }, // V
	/* T */ {false, false, false, false, false, false, false, true , true , true , true , true , true , true , true }, // T
	/* W */ {true , false, true , false, true , false, true , true , true , true , true , true , true , true , true }, // W
	/* Y */ {false, true , true , false, false, true , true , true , true , true , true , true , true , true , true }, // Y
	/* H */ {true , true , true , false, true , true , true , true , true , true , true , true , true , true , true }, // H
	/* K */ {false, false, false, true , true , true , true , true , true , true , true , true , true , true , true }, // K
	/* D */ {true , false, true , true , true , true , true , true , true , true , true , true , true , true , true }, // D
	/* B */ {false, true , true , true , true , true , true , true , true , true , true , true , true , true , true }, // B
	/* N */ {true , true , true , true , true , true , true , true , true , true , true , true , true , true , true }  // N
	//         A      C      M      G      R      S      V      T      W      Y      H      K      D      B      N
};

uint iupac_2_uint (const char & c) {
	switch (c) {
		case 'A':
			return 1;
		case 'C':
			return 2;
		case 'M':
			return 3;
		case 'G':
			return 4;
		case 'R':
			return 5;
		case 'S':
			return 6;
		case 'V':
			return 7;
		case 'T':
			return 8;
		case 'W':
			return 9;
		case 'Y':
			return 10;
		case 'H':
			return 11;
		case 'K':
			return 12;
		case 'D':
			return 13;
		case 'B':
			return 14;
		case 'N':
			return 15;
		default:
			return 0;
	}
	return 0;
}

bool iupac_comp(const char & a, const char & b) {
	uint iupac_a = iupac_2_uint(a);
	uint iupac_b = iupac_2_uint(b);

	if (iupac_a == 0) {
		cerr << a << " is not a IUPAC symbol !" << endl;
		exit (1);
	} else if (iupac_b == 0) {
		cerr << b << " is not a IUPAC symbol !" << endl;
		exit (1);
	}

	return matrix[iupac_a -1][iupac_b -1];
}


// int main () {
// 	cout << iupac_comp('A', 'C') << endl;
// 	cout << iupac_comp('A', 'M') << endl;
// 	cout << iupac_comp('N', 'C') << endl;

// 	if (iupac_comp('N', 'C'))
// 		cout << "Yeah !" << endl;

// 	return 0;
// }
