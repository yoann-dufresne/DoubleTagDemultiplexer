#include <iostream>
#include <fstream>
#include <vector>

#include "demultiplexing.hpp"
#include "Parser.hpp"
#include "edit.hpp"


using namespace std;


bool mistag = false;
bool trim = false;
ofstream mistag_r1;
ofstream mistag_r2;


vector<int> find_primers (const vector<Sequence> & primers, const string & r1, const string & r2, const uint errors) {
	vector<int> locations;
	string reads[2]; reads[0] = r1; reads[1] = r2;

	for (string & read : reads) {
		int primer_idx = -1;
		uint min_dist = errors+1;

		int primer_selected = -1;
		int position = -1;

		for (auto & primer : primers) {
			primer_idx++;

			string prim_seq = primer.sequence;

			// Get the begining of the sequence to perform alignment
			string read_start = read.substr(0, prim_seq.length());
			// Align
			compute_matrix(prim_seq, read_start);

			uint dist = get_distance();
			if (dist <= errors && dist < min_dist) {
				// Add the primer found in the result.
				primer_selected = primer_idx;
				min_dist = dist;

				// Compute the position to trim
				string symbols = get_align_symbols(prim_seq, read_start);
				// Compute global gaps to trim at the right place
				int gaps = 0;
				for (uint idx=0 ; idx<symbols.length() ; idx++)
					if (symbols[idx] == '_')
						gaps--;
					else if (symbols[idx] == '-')
						gaps++;
				position = prim_seq.length()+gaps;
			}
		}

		locations.push_back(primer_selected);
		locations.push_back(position);
	}

	return locations;
}


void demux (string r1_filename, string r2_filename,
	map<string, Experiment> exps,
	vector<Sequence> primers, uint errors) {

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

		// Compute distances
		vector<int> found = find_primers(primers, r1.sequence, r2.sequence, errors);

		int idx1 = found[0];
		int idx2 = found[2];

		if (idx1 != -1 && idx2 != -1) {
			nbFound++;

			Sequence p1 = primers[idx1];
			Sequence p2 = primers[idx2];

			// Triming primers
			if (trim) {
				r1.sequence = r1.sequence.substr(found[1]);
				r1.quality = r1.quality.substr(found[1]);
				r2.sequence = r2.sequence.substr(found[3]);
				r2.quality = r2.quality.substr(found[3]);
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
			mistag_r1 << (idx1 == -1 ? "unknown" : primers[idx1].header) << endl;
			mistag_r1 << r1.sequence << endl << "+" << endl << r1.quality << endl;

			// R2 read
			mistag_r2 << r2.header << ";tag:";
			mistag_r2 << (idx2 == -1 ? "unknown" : primers[idx2].header) << endl;
			mistag_r2 << r2.sequence << endl << "+" << endl << r2.quality << endl;
		}

		total++;
	}

	// Close files
	for (auto it=exps.begin() ; it!=exps.end() ; it++)
		it->second.closeFile();

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
