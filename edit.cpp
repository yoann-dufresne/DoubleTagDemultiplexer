#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>

#include "IUPAC.hpp"
#include "edit.hpp"

#define MAX_SIZE 250

// Constants
vector<vector<uint8_t> > matrix(MAX_SIZE+1, vector<uint8_t>(MAX_SIZE+1, 0));
uint lines=0, columns=0;


string matrix_2_string () {
	stringstream ss;

	for (uint line=0 ; line<lines ; line++) {
		for (uint column=0 ; column<columns ; column++) {
			ss << "\t" << (uint)matrix[line][column];
		}
		ss << endl;
	}

	return ss.str();
}

void compute_matrix(const string & s1, const string & s2) {
	// Too long words
	if (s1.length() >= MAX_SIZE || s2.length() >= MAX_SIZE) {
		cerr << "Words too long. MAX_SIZE: " << MAX_SIZE << endl;
		exit(1);
	}

	// Update matrix values
	lines = s1.length()+1;
	columns = s2.length()+1;

	// init
	for (uint x=0 ; x<lines ; x++)
		matrix[x][0] = x;
	for (uint y=0 ; y<columns ; y++)
		matrix[0][y] = y;

	// Fill the matrix
	for (uint x=0 ; x<s1.length() ; x++) {
		for (uint y=0 ; y<s2.length() ; y++) {
			// compute values from insertions deletion
			int del = matrix[x][y+1] + 1;
			int ins = matrix[x+1][y] + 1;

			// Compute the diagonal value
			int diag = matrix[x][y];
			// Add 1 if mismatch
			diag += iupac_comp(s1[x], s2[y]) ? 0 : 1;

			matrix[x+1][y+1] = min(min(ins, del), diag);
		}
	}
}

uint get_distance () {
	return matrix[lines-1][columns-1];
}

string get_align_symbols (const string & s1, const string & s2) {
	uint row = s1.length();
	uint col = s2.length();
	string align = "";

	while (row != 0 || col != 0) {
		// base conditions
		if (row == 0) {
			align = "-" + align;
			col--;
		} else if (col == 0) {
			align = "_" + align;
			row--;
		} else {
			// Substitution
			if (!iupac_comp(s1[row-1], s2[col-1]) && matrix[row-1][col-1] == matrix[row][col]-1) {
				align = " " + align;
				row--; col--;
			}
			// Match
			else if (iupac_comp(s1[row-1], s2[col-1]) && matrix[row-1][col-1] == matrix[row][col]) {
				align = "|" + align;
				row--; col--;
			}
			// Insertion
			else if (matrix[row-1][col] == matrix[row][col]-1) {
				align = "_" + align;
				row--;
			}
			// Deletion
			else if (matrix[row][col-1] == matrix[row][col]-1) {
				align = "-" + align;
				col--;
			}
			// Oops !
			else {
				cerr << "Woops ! Something wrong appened in the alignment matrix !" << endl;
				cerr << "Here are the coordinates and the matrix:" << endl;
				cerr << "row: " << row << "   col: " << col << endl;
				cerr << matrix_2_string();
				exit(1);
			}
		}
	}

	return align;
}

string get_alignment (const string s1, const string s2) {
	string symbols = get_align_symbols(s1, s2);

	// Rewrited words
	string w1 = "", w2 = "";
	uint s1_idx = 0, s2_idx = 0;

	for (uint idx=0 ; idx<symbols.length() ; idx++) {
		// Sub / Match
		if (symbols[idx] == ' ' || symbols[idx] == '|') {
			w1 += s1[s1_idx++];
			w2 += s2[s2_idx++];
		}
		// Insertion
		else if (symbols[idx] == '_') {
			w1 += s1[s1_idx++];
			w2 += "-";
		}
		// Deletion
		else {
			w1 += "-";
			w2 += s2[s2_idx];
		}
	}

	return w1 + "\n" + symbols + "\n" + w2 + "\n";
}


// int main(int argc, char* argv[]) {
// 	string w1 = string(argv[1]);
// 	string w2 = string(argv[2]);


// 	compute_matrix (w1, w2);
// 	cout << w1 << ' ' << w2 << endl;
// 	cout << matrix_2_string();
// 	cout << "dist: " << get_distance() << endl;
// 	cout << get_alignment(w1, w2);

// 	return 0;
// }
