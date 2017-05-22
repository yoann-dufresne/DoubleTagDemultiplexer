#include <string>
#include <fstream>

#include "Sequence.hpp"

#ifndef EXP_H
#define EXP_H

class Experiment {
public:
	std::string run;
	std::string sample;
	std::string fwd_name;
	std::string rev_name;

	std::string outname_fwd;
	std::string outname_rev;
	std::ofstream out_fwd;
	std::ofstream out_rev;

	Experiment() {};
	Experiment(std::string run, std::string sample, std::string fwd_tag, std::string rev_tag,
		std::string out_fwd, std::string out_rev) {
		this->run = run;
		this->fwd_name = fwd_tag;
		this->rev_name = rev_tag;
		this->outname_fwd = out_fwd;
		this->outname_rev = out_rev;
	};

	Experiment (const Experiment & e) {
		this->run = e.run;
		this->sample = e.sample;
		this->fwd_name = e.fwd_name;
		this->rev_name = e.rev_name;
		this->outname_fwd = e.outname_fwd;
		this->outname_rev = e.outname_rev;
	}

	Experiment& operator=(const Experiment& e) {
		this->run = e.run;
		this->sample = e.sample;
		this->fwd_name = e.fwd_name;
		this->rev_name = e.rev_name;
		this->outname_fwd = e.outname_fwd;
		this->outname_rev = e.outname_rev;

		return *this;
	}

	void addReads (Sequence fwd, Sequence rev);
	void closeFile ();
};

#endif
