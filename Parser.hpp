#include <string>
#include <fstream>

#include "Sequence.hpp"

#ifndef PARSER_H
#define PARSER_H

class Parser {
	std::string nextLine;
	std::ifstream stream;

	Sequence onFileError();
	Sequence nextFasta();
	Sequence nextFastq();

public:
	std::string filename;
	bool verbose;

	/* Init a fasta/q parser on a file */
	Parser (std::string filename);

	/* True if there is another sequence in the file */
	bool hasNext ();

	/* Read the next sequence fatsa/q)
	 * The quality will be empty if it's a fasta file
	 * Return an empty sequence if the file have problem or when eof is encountered
	 */
	Sequence nextSequence ();
};

#endif
