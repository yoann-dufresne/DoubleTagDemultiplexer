#include <iostream>
#include <fstream>
#include <vector>
#include <stdexcept>

#include "demultiplexing.hpp"
#include "Parser.hpp"
#include "edit.hpp"
#include "IUPAC.hpp"


using namespace std;


bool mistag = false;
bool trim = false;
ofstream mistag_r1;
ofstream mistag_r2;
ofstream empty;


int find_with_error (string & prim_seq, string & read_start, uint errors);
int find_0_error (string prim_seq, string read_start);


uint min_dist;
vector<int> find_primers (const vector<Sequence> & primers, const string & r1, const string & r2, const uint errors) {
	vector<int> locations;
	string reads[2]; reads[0] = r1; reads[1] = r2;

	for (string & read : reads) {
		int primer_idx = -1;
		min_dist = errors+1;

		int primer_selected = -1;
		int position = -1;

		for (auto & primer : primers) {
			primer_idx++;

			string prim_seq = primer.sequence;

			// Print an error on length problems
			if (read.length() < prim_seq.length()) {
				cerr << "Sequence too short regarding the primer" << endl;
				cerr << "Read:   " << read << endl;
				cerr << "Primer: " << prim_seq << endl;
				cerr << "Skipping..." << endl;
				continue;
			}

			// Get the begining of the sequence to perform alignment
			string read_start = read.substr(0, prim_seq.length());
			
			int tmp_position = -1;
			if (errors == 0)
				// strict align
				try {
					tmp_position = find_0_error(prim_seq, read_start);
				} catch (const invalid_argument & ia) {
					cerr << "Read:   " << read << endl;
					cerr << "Primer: " << prim_seq << endl;
					cerr << "Skipping..." << endl;
					continue;
				}
			else
				// Align with edit distance
				tmp_position = find_with_error (prim_seq, read_start, errors);

			if (tmp_position != -1) {
				position = tmp_position;
				primer_selected = primer_idx;
			}
		}

		locations.push_back(primer_selected);
		locations.push_back(position);
	}

	return locations;
}

int find_0_error (string prim_seq, string read_start) {
	for (uint idx=0 ; idx<prim_seq.length() ; idx++) {
		try {
			if (prim_seq[idx] != read_start[idx] && !iupac_comp(prim_seq[idx], read_start[idx]))
				return -1;
		} catch (const invalid_argument& ia) {
			cerr << "Wrong IUPAC symbol at idx " << idx << endl;
			throw ia;
		}
	}
	return prim_seq.length();
}

int find_with_error (string & prim_seq, string & read_start, uint errors) {
	// Align
	compute_matrix(prim_seq, read_start);

	uint dist = get_distance();
	if (dist <= errors && dist < min_dist) {
		// Add the primer found in the result.
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
		
		int position = prim_seq.length()+gaps;

		return position;
	}
	return -1;
}


void write_empty (Sequence & s1, Sequence & s2) {
	empty << ">" << s1.header << " from_R1" << endl;
	empty << s1.sequence << endl;
	empty << ">" << s2.header << " from_R2" << endl;
	empty << s2.sequence << endl;
}


void write_mistag (Sequence & s1, Sequence & s2, string tag1, string tag2) {
	// R1 read
	mistag_r1 << "@" << s1.header << ";tag:" << tag1 << endl;
	mistag_r1 << s1.sequence << endl << "+" << endl << s1.quality << endl;

	// R2 read
	mistag_r2 << "@" << s2.header << ";tag:" << tag2 << endl;
	mistag_r2 << s2.sequence << endl << "+" << endl << s2.quality << endl;
}


void demux (string r1_filename, string r2_filename,
	map<string, Experiment> exps,
	vector<Sequence> primers, uint errors, bool check_end, uint min_length) {

	// For all the paired end reads
	Parser r1_parse (r1_filename), r2_parse(r2_filename);
	r1_parse.verbose = true;
	r2_parse.verbose = true;
	int nbFound = 0;
	int exp_count = 0;
	int empty_count = 0;
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

				if (check_end) {
					// Create reverse complement of each sequence
					Sequence * r1rc = r1.revcomp();
					Sequence * r2rc = r2.revcomp();
					// Check for the primers at the beginning of the reverse complemented read
					found = find_primers(primers, r1rc->sequence, r2rc->sequence, errors);
					free(r1rc);
					free(r2rc);

					// Trim the reverse complement primer if founded.
					if (found[0] != -1) {
						r1.sequence = r1.sequence.substr(0, r1.sequence.length()-found[1]-1);
						if (r1.sequence.length() == 0)
							cerr << "Problem with empty read " << r1.header << endl;
						r1.quality = r1.quality.substr(0, r1.quality.length()-found[1]-1);
					}
					if (found[2] != -1) {
						r2.sequence = r2.sequence.substr(0, r2.sequence.length()-found[3]-1);
						if (r2.sequence.length() == 0)
							cerr << "Problem with empty read " << r2.header << endl;
						r2.quality = r2.quality.substr(0, r2.quality.length()-found[3]-1);
					}
				}
			}

			// Extract empty reads
			if (mistag && (r1.sequence.length() < min_length || r2.sequence.length() < min_length)) {
				write_empty(r1, r1);
				empty_count++;
			} else {
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
					} else {
						write_mistag(r1, r2, primers[idx1].header, primers[idx2].header);
					}
				}
			}
		} else if (mistag) {
			// Save the mistag
			write_mistag(
				r1, r2, 
				idx1 == -1 ? "unknown" : primers[idx1].header,
				idx2 == -1 ? "unknown" : primers[idx2].header
			);
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
		empty.close();
	}

	cout << "Input reads: " << total << endl;
	cout << "Empty reads: " << empty_count << endl;
	cout << "No primer found: " << (total-nbFound-empty_count) << endl;
	cout << "Unasignable: " << (total-exp_count) << endl;
}

void activate_mistags (string out_dir, string run_name) {
	mistag = true;

	// Add the directory mark
	if (out_dir[out_dir.length() - 1] != '/')
		out_dir += "/";

	// Open mistag files
	mistag_r1.open(out_dir + run_name + "_mistag_R1.fastq");
	mistag_r2.open(out_dir + run_name + "_mistag_R2.fastq");
	empty.open(out_dir + run_name + "_empty.fasta");
}

void activate_triming () {
	trim = true;
}
