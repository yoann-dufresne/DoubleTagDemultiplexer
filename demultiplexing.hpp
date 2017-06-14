#include <string>
#include <map>

#include "Experiment.hpp"
#include "Sequence.hpp"

void demux (std::string r1_filename, std::string r2_filename,
	std::map<std::string, Experiment> exps, std::map<std::string, Sequence> primers);

void activate_mistags (std::string out_dir);
void activate_triming ();
