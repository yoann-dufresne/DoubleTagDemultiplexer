#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>
#include <map>
#include <cstdio>

#include <sys/stat.h>

#include <boost/filesystem.hpp>

#include "Sequence.hpp"
#include "Parser.hpp"
#include "Experiment.hpp"

#include "demultiplexing.hpp"

using namespace std;
namespace fs=boost::filesystem;

void printHelp () {
	cout << "./dtd -r1 r1_filename -r2 r2_filename \
-o oligos_filename -e experiments [-d output_directory]" << endl;
	cout << "Command line description:" << endl;
	cout << "\t\e[1m-r1\e[0m: The FASTQ file containing all the reads R1 from \
a paired end amplicon sequencing." << endl;
	cout << "\t\e[1m-r2\e[0m: The FASTQ file containing all the reads R2 from \
a paired end amplicon sequencing." << endl;
	cout << "\t\e[1m-o or -oligos\e[0m: The FASTA file containing all the tag-primer pairs. The \
header line must contain the name of the sequence using the convention \
<primer_name>-<tag_name>. For example F1-C. The primer can contain the '-' \
character but nor the tag name." << endl;
	// cout << "\t\e[1m-t or -tag-size\e[0m: The size of the tags in bp." << endl;
	cout << "\t\e[1m-e or -experiments\e[0m: The csv file containing the tags by experiment. \
The csv must have at least 4 columns named run, sample, forward and reverse. \
These columns corresponds respectively to the run name, the sample name, the forward \
taged primer name of the experiment and the reverse." << endl;
	cout << "\t\e[1m-d or -directory\e[0m: The directory where all the output FASTQ by \
by experiment will be created. By default this is the current directory." << endl;
	cout << "/t/e[1m-m or -mistag\e[0m: Two files named mistag_R1.fastq and \
mistag_R2.fastq will be created with non assigned read pairs. The primer names will be \
added to the end of the read headers." << endl;
}

map<string, Experiment> parse_experiments(string exp_filename, string out_dir);
map<string, Sequence> parseTaggedPrimers (string filename);

int main (int argc, char *argv[]) {
	string r1_filename;
	string r2_filename;
	string exp_filename;
	string oligos_filename;
	int tag_size = -1;
	string out_dir = "./";
	bool mistags = false;

	if (argc == 1) {
		printHelp();
		return 0;
	}

	/* --- Arguments parsing --- */
	for (int idx=1 ; idx<argc ; idx++) {
		string arg (argv[idx]);

		if (arg == "-h" || arg == "-help") {
			printHelp();
			return 0;
		} else if (arg == "-r1") {
			arg = string(argv[++idx]);
			r1_filename = arg;
		} else if (arg == "-r2") {
			arg = string(argv[++idx]);
			r2_filename = arg;
		} else if (arg == "-o" || arg == "-oligos") {
			arg = string(argv[++idx]);
			oligos_filename = arg;
		} else if (arg == "-t" || arg == "-tag-size") {
			arg = string(argv[++idx]);
			istringstream buffer(arg);
			buffer >> tag_size;
		} else if (arg == "-e" || arg == "-experiment") {
			arg = string(argv[++idx]);
			exp_filename = arg;
		} else if (arg == "-d" || arg == "-directory") {
			arg = string(argv[++idx]);
			out_dir = arg;
			if (out_dir[0] != '/')
				out_dir = "./" + out_dir;
		} else if (arg == "-m" || arg == "-mistag") {
			mistags = true;
		} else {
			cerr << "No argument called " << arg << endl;
			return 0;
		}
	}

	if (r1_filename == "" || r2_filename == "" || exp_filename == "" || oligos_filename == "") {
		cerr << "Wrong command line" << endl;
		return 1;
	}

	/* --- Loading the tags --- */
	auto exps = parse_experiments(exp_filename, out_dir);
	auto oligos = parseTaggedPrimers(oligos_filename);

	/* --- Demultiplexing --- */
	// Create directory
	fs::path dir = fs::path(out_dir);
	if (! fs::exists(dir))
		if (mkdir(out_dir.c_str(), 0777)) {
			cerr << "Impossible to create directory " << out_dir << endl;
			exit(2);
		}

	// Mistage saving
	if (mistags)
		activate_mistags (out_dir);
	// Demultiplex
	demux (r1_filename, r2_filename, exps, oligos);
	for (auto it=exps.begin() ; it!=exps.end() ; it++) {
		it->second.closeFile();
	}

	return 0;
}


vector<string> split (string word, string delim) {
	vector<string> words;

	size_t pos = 0;
	string token;
	while ((pos = word.find(delim)) != string::npos) {
		token = word.substr(0, pos);
		words.push_back(token);
		word.erase(0, pos + delim.length());
	}
	words.push_back(word);

	return words;
}

map<char, string> getIupac () {
	map<char, string> iupac;

	iupac['A'] = "A";
	iupac['C'] = "C";
	iupac['T'] = "T";
	iupac['G'] = "G";
	iupac['R'] = "AG";
	iupac['Y'] = "CT";
	iupac['S'] = "GC";
	iupac['W'] = "AT";
	iupac['K'] = "GT";
	iupac['M'] = "AC";
	iupac['B'] = "CGT";
	iupac['D'] = "AGT";
	iupac['H'] = "ACT";
	iupac['V'] = "ACG";
	iupac['N'] = "ACGT";

	return iupac;
}

/* Recursively split the sequence into multiple sequences, each time that a degenerate
 * nucleotide is encountered
 */
vector<Sequence> expandWithIupacCode (Sequence & seq, map<char, string> & iupac, int startIdx) {
	vector<Sequence> expanded;

	for (uint idx=startIdx ; idx<seq.sequence.size() ; idx++) {
		auto iVal = iupac[seq.sequence[idx]];

		if (iVal.size() > 1) {
			for (char c : iVal) {
				Sequence copy = Sequence(seq);
				copy.sequence[idx] = c;
				vector<Sequence> recur = expandWithIupacCode(copy, iupac, idx+1);
				expanded.insert(expanded.end(), recur.begin(), recur.end());
			}
		}
	}

	if (expanded.size() == 0)
		expanded.push_back(seq);

	return expanded;
}

map<string, Experiment> parse_experiments(string exp_filename, string out_dir) {
	ifstream file (exp_filename);

	map<string, Experiment> exps;

	// Wrong file of path
	if (!file.is_open()) {
		cerr << "Impossible to open " << exp_filename << endl;
		exit(1);
	}

	// Read the header line
	string line;
	getline(file, line);
	vector<string> header = split(line, ",");

	int run_idx = -1;
	int sample_idx = -1;
	int r1_idx = -1;
	int r2_idx = -1;

	for (uint idx=0 ; idx<header.size() ; idx++) {
		string word = header[idx];

		if (word == "run")
			run_idx = idx;
		else if (word == "sample")
			sample_idx = idx;
		else if (word == "forward")
			r1_idx = idx;
		else if (word == "reverse")
			r2_idx = idx;
	}

	// If the csv is not correctly formated
	if (run_idx == -1)
		cerr << "No field run in " << exp_filename << endl;
	if (sample_idx == -1)
		cerr << "No field sample in " << exp_filename << endl;
	if (r1_idx == -1)
		cerr << "No field forward in " << exp_filename << endl;
	if (r2_idx == -1)
		cerr << "No field reverse in " << exp_filename << endl;
	if(run_idx == -1 || sample_idx == -1 || r1_idx == -1 || r2_idx == -1)
		exit(1);

	while (file.peek() != EOF) {
		getline(file, line);
		vector<string> values = split(line, ",");

		if (values.size() < 4)
			continue;

		string base_filename = out_dir + (out_dir[out_dir.size()-1] == '/' ? "" : "/")
			+ values[run_idx] + "_" + values[sample_idx];
		Experiment e(values[run_idx], values[sample_idx], values[r1_idx], values[r2_idx],
			base_filename + "_fwd.fastq", base_filename + "_rev.fastq");
		exps[e.fwd_name + e.rev_name] = e;
	}

	return exps;
}


map<string, Sequence> parseTaggedPrimers (string filename) {
	map<string, Sequence> oligos;

	uint minSize = 1000000000;
	vector<Sequence> seqs;
	Parser fasta (filename);
	while (fasta.hasNext()) {
		Sequence seq = fasta.nextSequence();
		if (seq.sequence.size() < minSize)
			minSize = seq.sequence.size();
		seqs.push_back(seq);
	}

	auto iupac = getIupac();
	for (Sequence seq : seqs) {
		vector<Sequence> expanded = expandWithIupacCode(seq, iupac, 0);

		for (Sequence & exp : expanded) {
			oligos[exp.sequence.substr(0, minSize)] = seq;
		}
	}

	return oligos;
}

