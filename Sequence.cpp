#include "Sequence.hpp"

using namespace std;

Sequence::Sequence () {}

Sequence::Sequence (const Sequence& seq) {
	this->header = seq.header;
	this->sequence = seq.sequence;
	this->quality = seq.quality;
}
