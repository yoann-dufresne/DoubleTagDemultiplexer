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
-p primers_filename -e experiments [-d output_directory]" << endl;
	cout << "Command line description:" << endl;
	cout << "\t\e[1m-r1\e[0m: The FASTQ file containing all the reads R1 from \
a paired end amplicon sequencing." << endl;
	cout << "\t\e[1m-r2\e[0m: The FASTQ file containing all the reads R2 from \
a paired end amplicon sequencing." << endl;
	cout << "\t\e[1m-p or -primers\e[0m: The FASTA file containing all the tag-primer pairs. The \
header line must contain the name of the sequence using the convention \
<primer_name>-<tag_name>. For example F1-C. The primer can contain the '-' \
character but nor the tag name." << endl;
	// cout << "\t\e[1m-t or -tag-size\e[0m: The size of the tags in bp." << endl;
	cout << "\t\e[1m-l or -libraries\e[0m: The csv file containing the tags for the libraries and samples. \
The csv must have at least 4 columns named run, sample, forward and reverse. \
These columns corresponds respectively to the run name, the sample name, the forward \
taged primer name of the experiment and the reverse." << endl;
	cout << "\t\e[1m-d or -directory\e[0m: The directory where all the output FASTQ by \
by experiment will be created. By default this is the current directory." << endl;
	cout << "\t\e[1m-m or -mistag\e[0m: Two files named mistag_R1.fastq and \
mistag_R2.fastq will be created with non assigned read pairs. The primer names will be \
added to the end of the read headers." << endl;
	cout << "\t\e[1m-t or -trim\e[0m: Trim the primers from the sequences." << endl;
	cout << "\t\e[1m-e or -errors\e[0m: A value representing the number of errors allowed in the primers." << endl;
	cout << "\t\e[1m-rl or -restrict-library\e[0m: Restrict the demultiplexing to only one library." << endl;
}

map<string, Experiment> parse_experiments(string exp_filename, string out_dir, string restriction);
vector<Sequence> parseTaggedPrimers (string filename);

int main (int argc, char *argv[]) {
	string r1_filename;
	string r2_filename;
	string exp_filename;
	string primers_filename;
	string out_dir = "./";
	string restricted = "";
	bool mistags = false;
	bool trim = false;
	uint e = 0;

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
		} else if (arg == "-p" || arg == "-primers") {
			arg = string(argv[++idx]);
			primers_filename = arg;
		} else if (arg == "-l" || arg == "-libraries") {
			arg = string(argv[++idx]);
			exp_filename = arg;
		} else if (arg == "-d" || arg == "-directory") {
			arg = string(argv[++idx]);
			out_dir = arg;
			if (out_dir[0] != '/')
				out_dir = "./" + out_dir;
		} else if (arg == "-m" || arg == "-mistag") {
			mistags = true;
		} else if (arg == "-t" || arg == "-trim") {
			trim = true;
		} else if (arg == "-e" || arg == "-errors") {
			e = atoi(argv[++idx]);
		} else if (arg == "-rl" || arg == "-restrict-library") {
			restricted = string(argv[++idx]);
		} else {
			cerr << "No argument called " << arg << endl;
			return 1;
		}
	}

	if (r1_filename == "" || r2_filename == "" || exp_filename == "" || primers_filename == "") {
		cerr << "Wrong command line" << endl;
		return 1;
	}

	/* --- Loading the tags --- */
	auto exps = parse_experiments(exp_filename, out_dir, restricted);
	auto primers = parseTaggedPrimers(primers_filename);

	/* --- Demultiplexing --- */
	// Create directory
	fs::path dir = fs::path(out_dir);
	if (! fs::exists(dir))
		if (mkdir(out_dir.c_str(), 0777)) {
			cerr << "Impossible to create directory " << out_dir << endl;
			exit(2);
		}

	// Flags activation
	if (mistags) {
		if (restricted != "")
			activate_mistags (out_dir, restricted);
		else
			activate_mistags (out_dir, exps.begin()->second.run);
	}
	if (trim)
		activate_triming ();

	// Demultiplex
	demux (r1_filename, r2_filename, exps, primers, e);

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

map<string, Experiment> parse_experiments(string exp_filename, string out_dir, string restriction) {
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

		if (values.size() < 4 || (restriction != "" && values[run_idx] != restriction))
			continue;

		string short_filename = exp_filename.substr(0, exp_filename.find_last_of("."));
		if (short_filename.find_last_of("/") != string::npos)
			short_filename = short_filename.substr(short_filename.find_last_of("/"));

		string base_filename =  out_dir + (out_dir[out_dir.size()-1] == '/' ? "" : "/")
			+ short_filename + '_' + values[run_idx] + "_" + values[sample_idx];
		Experiment e(values[run_idx], values[sample_idx], values[r1_idx], values[r2_idx],
			base_filename + "_fwd.fastq", base_filename + "_rev.fastq");
		
		// Error if primer pair is already used
		if (exps.find(e.fwd_name + e.rev_name) != exps.end()) {
			cerr << "The primer pair (" << e.fwd_name << " / " << e.rev_name << ") is used a multiple time in the same run." << endl;
			cerr << "Due to this fact, experiments are abiguous and the demultiplexer can't sort reads." << endl;
			cerr << "Please verify your tag to samples file." << endl;
			exit(1);
		}

		exps[e.fwd_name + e.rev_name] = e;
	}

	return exps;
}


vector<Sequence> parseTaggedPrimers (string filename) {
	vector<Sequence> primers;

	Parser fasta (filename);
	while (fasta.hasNext()) {
		Sequence seq = fasta.nextSequence();
		primers.push_back(seq);
	}

	return primers;
}

