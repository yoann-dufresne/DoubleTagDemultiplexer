#include <iostream>
#include <algorithm>

#include "Parser.hpp"

using namespace std;


Parser::Parser (string filename) {
	this->filename = filename;
	this->stream.open(filename);
	if (!this->stream.is_open()) {
		cerr << "Impossible to load " << filename << endl;
		exit(1);
	}
	this->verbose = false;
	getline(this->stream, this->nextLine);
}

Sequence Parser::nextSequence () {
	// File ended
	if (! this->stream.is_open())
		return Sequence();

	if (this->nextLine[0] == '>')
		return this->nextFasta();
	else if (this->nextLine[0] == '@')
		return this->nextFastq();
	else {
		return this->onFileError();
	}
}

Sequence Parser::onFileError () {
	cerr << "Unexpected line:" << endl;
	cerr << this->nextLine << endl;
	cerr << "Skipping..." << endl;

	this->stream.close();
	return Sequence();
}

bool Parser::hasNext () {
	return this->stream.is_open();
}

Sequence Parser::nextFasta () {
	Sequence seq;

	// Header
	this->nextLine.erase(this->nextLine.begin());
	seq.header = this->nextLine;
	// Read next line
	if (!getline(this->stream, this->nextLine)) {
		this->stream.close();
		return Sequence();
	}

	// Sequence
	while (this->stream.is_open() && this->nextLine[0] != '>') {
		seq.sequence += this->nextLine;
		seq.sequence.erase(remove(seq.sequence.begin(), seq.sequence.end(), '\n'), seq.sequence.end());

		if (!getline(this->stream, this->nextLine)) {
			this->stream.close();
			return seq;
		}
	}

	return seq;
}

Sequence Parser::nextFastq () {
	Sequence seq;

	// Header
	this->nextLine.erase(this->nextLine.begin());
	seq.header = this->nextLine;
	// Read next line
	if (!getline(this->stream, this->nextLine)) {
		this->stream.close();
		return Sequence();
	}

	// Sequence
	while (this->nextLine[0] != '+') {
		seq.sequence += this->nextLine;

		if (!getline(this->stream, this->nextLine))
			return this->onFileError();
	}

	// Read the line after the '+'
	if (!getline(this->stream, this->nextLine))
		return this->onFileError();

	// Sequence
	seq.quality += this->nextLine;

	if (!getline(this->stream, this->nextLine)) {
		this->stream.close();
		return seq;
	}

	return seq;
}
