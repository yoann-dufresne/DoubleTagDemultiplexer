#include <string>


#ifndef SEQ_H
#define SEQ_H

class Sequence {
public:
	std::string header;
	std::string sequence;
	std::string quality;

	Sequence();
	Sequence(const Sequence& seq);

	Sequence * revcomp();
};


#endif
