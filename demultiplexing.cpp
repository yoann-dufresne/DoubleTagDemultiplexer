#include <iostream>
#include <fstream>

#include "demultiplexing.hpp"
#include "Parser.hpp"


using namespace std;


bool mistag = false;
bool trim = false;
ofstream mistag_r1;
ofstream mistag_r2;


void demux (string r1_filename, string r2_filename,
	map<string, Experiment> exps, map<string, Sequence> primers) {

	int split_size = primers.begin()->first.size();

	// For all the paired end reads
	Parser r1_parse (r1_filename), r2_parse(r2_filename);
	r1_parse.verbose = true;
	r2_parse.verbose = true;
	int nbFound = 0;
	int exp_count = 0;
	int total = 0;
	while (r1_parse.hasNext() && r2_parse.hasNext()) {
		Sequence r1, r2;
		r1 = r1_parse.nextSequence();
		r2 = r2_parse.nextSequence();

		// Primers
		string s1, s2;
		s1 = r1.sequence.substr(0, split_size);
		s2 = r2.sequence.substr(0, split_size);

		auto it1 = primers.find(s1);
		auto it2 = primers.find(s2);
		if (it1 != primers.end() && it2 != primers.end()) {
			nbFound++;

			Sequence p1 = it1->second;
			Sequence p2 = it2->second;

			// Triming primers
			if (trim) {
				r1.sequence = r1.sequence.substr(p1.sequence.length());
				r2.sequence = r2.sequence.substr(p2.sequence.length());
			}

			// R1 == fwd
			auto it_fwd = exps.find(p1.header+p2.header);
			if (it_fwd != exps.end()) {
				exp_count++;
				it_fwd->second.addReads(r1, r2);
			} else {
				auto it_rev = exps.find(p2.header+p1.header);
				if (it_rev != exps.end()) {
					exp_count++;
					it_rev->second.addReads(r2, r1);
				}
			}			
		} else if (mistag) {
			// Save the mistag

			// R1 read
			mistag_r1 << r1.header << ";tag:";
			mistag_r1 << (it1 == primers.end() ? "unknown" : it1->second.header) << endl;
			mistag_r1 << r1.sequence << endl << "+" << endl << r1.quality << endl;

			// R2 read
			mistag_r2 << r2.header << ";tag:";
			mistag_r2 << (it2 == primers.end() ? "unknown" : it2->second.header) << endl;
			mistag_r2 << r2.sequence << endl << "+" << endl << r2.quality << endl;
		}

		total++;
	}

	// Close mistag files
	if (mistag) {
		mistag_r1.close();
		mistag_r2.close();
	}

	cout << "Input reads: " << total << endl;
	cout << "Primers found: " << nbFound << endl;
	cout << "Unasignable: " << (total-exp_count) << endl;
}

void activate_mistags (string out_dir) {
	mistag = true;

	// Add the directory mark
	if (out_dir[out_dir.length() - 1] != '/')
		out_dir += "/";

	// Open mistag files
	mistag_r1.open(out_dir + "mistag_R1.fastq");
	mistag_r2.open(out_dir + "mistag_R2.fastq");
}

void activate_triming () {
	trim = true;
}
