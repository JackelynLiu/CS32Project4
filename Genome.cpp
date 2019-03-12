#include "provided.h"
#include <string>
#include <vector>
#include <iostream>
#include <istream>
#include <fstream>
using namespace std;

class GenomeImpl
{
public:
	GenomeImpl(const string& nm, const string& sequence);
	static bool load(istream& genomeSource, vector<Genome>& genomes);
	int length() const;
	string name() const;
	bool extract(int position, int length, string& fragment) const;
private:
	string m_name, m_sequence;
};

GenomeImpl::GenomeImpl(const string& nm, const string& sequence)
	:m_name(nm), m_sequence(sequence)
{}

bool GenomeImpl::load(istream& genomeSource, vector<Genome>& genomes)
{
	string name, sequence, line;
	while (getline(genomeSource, line))
	{
		if (line[0] == '>')
		{
			if (name.size() != 0 && sequence.size() != 0)
			{
				genomes.push_back(Genome(name, sequence));
				name = "";
				sequence = "";
			}
			if (line[1] == ' ' || line[1] == '\n') return false;
			for (int i=1;i != '\n';i++)
				name += line[i];
		}
		else
		{
			for (int i = 0; i != line.size(); i++)
			{
				char n = line[i];
				switch (n)
				{
				case 'A':
				case 'C':
				case 'T':
				case 'G':
				case 'N':
					sequence += line[i];
					break;
				default:
					return false;
				}
			}
		}
	}

	if (name.size() != 0 && sequence.size() != 0)
	{
		genomes.push_back(Genome(name, sequence));
		name, sequence = "";
	}
	return true;
}

int GenomeImpl::length() const
{
	return m_sequence.length();
}

string GenomeImpl::name() const
{
	return m_name;
}

bool GenomeImpl::extract(int position, int length, string& fragment) const
{
	if (m_sequence.size() - position < length) return false;
	else
	{
		fragment = m_sequence.substr(position, length);
		return true;
	}
}

//******************** Genome functions ************************************

// These functions simply delegate to GenomeImpl's functions.
// You probably don't want to change any of this code.

Genome::Genome(const string& nm, const string& sequence)
{
	m_impl = new GenomeImpl(nm, sequence);
}

Genome::~Genome()
{
	delete m_impl;
}

Genome::Genome(const Genome& other)
{
	m_impl = new GenomeImpl(*other.m_impl);
}

Genome& Genome::operator=(const Genome& rhs)
{
	GenomeImpl* newImpl = new GenomeImpl(*rhs.m_impl);
	delete m_impl;
	m_impl = newImpl;
	return *this;
}

bool Genome::load(istream& genomeSource, vector<Genome>& genomes)
{
	return GenomeImpl::load(genomeSource, genomes);
}

int Genome::length() const
{
	return m_impl->length();
}

string Genome::name() const
{
	return m_impl->name();
}

bool Genome::extract(int position, int length, string& fragment) const
{
	return m_impl->extract(position, length, fragment);
}