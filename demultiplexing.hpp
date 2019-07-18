#include <string>
#include <map>
#include <vector>

#include "Experiment.hpp"
#include "Sequence.hpp"

void demux (std::string r1_filename, std::string r2_filename,
	std::map<std::string, Experiment> exps,
	std::vector<Sequence> primers, uint errors, bool check_end, uint min_length);

void open_empty_file(std::string out_dir, std::string run_name);
void activate_mistags (std::string out_dir, std::string run_name);
void activate_triming ();
