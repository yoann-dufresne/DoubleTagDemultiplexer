#include "Experiment.hpp"
#include <iostream>

using namespace std;

void Experiment::addReads (Sequence fwd, Sequence rev) {
	if (! this->out_fwd.is_open())
		this->out_fwd.open(this->outname_fwd);
	if (! this->out_rev.is_open())
		this->out_rev.open(this->outname_rev);

	// Save fwd read
	this->out_fwd << '>' << fwd.header << endl;
	this->out_fwd << fwd.sequence << endl;
	this->out_fwd << '+' << endl;
	this->out_fwd << fwd.quality << endl;

	// Save rev read
	this->out_rev << '>' << rev.header << endl;
	this->out_rev << rev.sequence << endl;
	this->out_rev << '+' << endl;
	this->out_rev << rev.quality << endl;
};

void Experiment::closeFile () {
	if (this->out_fwd.is_open())
		this->out_fwd.close();
	if (this->out_rev.is_open())
		this->out_rev.close();
};