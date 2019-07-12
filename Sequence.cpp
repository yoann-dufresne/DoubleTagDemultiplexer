#include "Sequence.hpp"
#include <sstream>

using namespace std;

Sequence::Sequence () {}

Sequence::Sequence (const Sequence& seq) {
	this->header = seq.header;
	this->sequence = seq.sequence;
	this->quality = seq.quality;
}

Sequence * Sequence::revcomp () {
	Sequence * s = new Sequence();
	// this->quality.copy(s->quality, sizeof this->quality);
	// reverse(this->quality.begin(), this->quality.end());

	stringstream ss;
	auto it=this->sequence.end();
	do {
		it--;
		char c = *it;
		switch(c) {
			case 'A':
			case 'a':
				ss << 'T';
				break;
			case 'C':
			case 'c':
				ss << 'G';
				break;
			case 'G':
			case 'g':
				ss << 'C';
				break;
			case 'T':
			case 't':
				ss << 'A';
				break;
			default:
				ss << 'N';
		}
	} while(it != sequence.begin());
	s->sequence = ss.str();

	return s;
}
