#include <string>
#include <fstream>

#ifndef EXP_H
#define EXP_H

class Experiment {
public:
	std::string run;
	std::string sample;
	std::string fwd_name;
	std::string rev_name;
	std::ofstream out_fwd;
	std::ofstream out_rev;

	Experiment() {};
	Experiment(std::string run, std::string sample, std::string fwd, std::string rev,
		std::string out_fwd, std::string out_rev) {
		this->run = run;
		this->fwd_name = fwd;
		this->rev_name = rev;
		this->out_fwd.open(out_fwd);
		this->out_rev.open(out_rev);
	};

	Experiment (const Experiment & e) {
		this->run = e.run;
		this->sample = e.sample;
		this->fwd_name = e.fwd_name;
		this->rev_name = e.rev_name;
		this->out_fwd.open(fwd_name);
		this->out_rev.open(rev_name);
	}

	Experiment& operator=(const Experiment& e) {
		this->run = e.run;
		this->sample = e.sample;
		this->fwd_name = e.fwd_name;
		this->rev_name = e.rev_name;
		this->out_fwd.open(fwd_name);
		this->out_rev.open(rev_name);

		return *this;
	}
};

#endif
