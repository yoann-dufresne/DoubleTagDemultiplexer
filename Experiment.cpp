#include "Experiment.hpp"
#include <iostream>

using namespace std;

void Experiment::addReads (Sequence fwd, Sequence rev) {
	if (! this->out_fwd.is_open())
		this->out_fwd.open(this->outname_fwd);
	if (! this->out_rev.is_open())
		this->out_rev.open(this->outname_rev);

	this->nbReads += 1;

	// Save fwd read
	this->out_fwd << (fwd.quality == "" ? '>' : '@') << fwd.header << endl;
	this->out_fwd << fwd.sequence << endl;
	if (fwd.quality != "") {
		this->out_fwd << '+' << endl;
		this->out_fwd << fwd.quality << endl;
	}

	// Save rev read
	this->out_rev << (rev.quality == "" ? '>' : '@') << rev.header << endl;
	this->out_rev << rev.sequence << endl;
	if (rev.quality != "") {
		this->out_rev << '+' << endl;
		this->out_rev << rev.quality << endl;
	}
};

void Experiment::closeFile () {
	if (this->nbReads == 0) {
		this->out_fwd.open(this->outname_fwd);
		this->out_rev.open(this->outname_rev);
	}

	// Print reads by sample
	cout << this->outname_fwd << " and " << this->outname_rev << ": " << this->nbReads << " reads" << endl;
	

	if (this->out_fwd.is_open())
		this->out_fwd.close();
	if (this->out_rev.is_open())
		this->out_rev.close();
};