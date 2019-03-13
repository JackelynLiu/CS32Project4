#include "provided.h"
#include "Trie.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
using namespace std;

class GenomeMatcherImpl
{
public:
	GenomeMatcherImpl(int minSearchLength);
	void addGenome(const Genome& genome);
	int minimumSearchLength() const;
	bool findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const;
	bool findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const;
private:
	struct Loc2pos
	{
		Loc2pos(int l, int p) :m_loc(l), m_pos(p) {}
		int m_loc, m_pos;
	};

	int m_minsearchlength;
	vector<Genome> m_genomes;
	Trie<Loc2pos> m_trieofDNA;
	int determineifExactMatchorSNiP(string first, string second) const;
	//void findGenomesWithThisDNAhelp(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const;
};

GenomeMatcherImpl::GenomeMatcherImpl(int minSearchLength)
	:m_minsearchlength(minSearchLength)
{}

void GenomeMatcherImpl::addGenome(const Genome& genome)
{
	m_genomes.push_back(genome);
	for (int i = 0; i < genome.length() - m_minsearchlength + 1; i++)
	{
		string s;
		if(genome.extract(i, m_minsearchlength, s))
			m_trieofDNA.insert(s, Loc2pos(m_genomes.size()-1, i));
	}
}

int GenomeMatcherImpl::minimumSearchLength() const
{
	return m_minsearchlength;
}

int GenomeMatcherImpl::determineifExactMatchorSNiP(string first, string second) const
{
	if (first.compare(second) == 0) return 1;		//returns 1 if exact
	bool mismatch_found = false;
	if (first[0] != second[0]) return 0;			//returns 0 if neither
	for (int i = 0; i < first.size(); i++)
	{
		if (first[i] != second[i])
		{
			if (mismatch_found)
				return 0;
			else mismatch_found = true;
		}
	}
	return 2;						//returns 2 if SNiP
}

bool GenomeMatcherImpl::findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const
{
	if (fragment.length() < minimumLength) return false;
	if (minimumLength < m_minsearchlength) return false;

	string tosearch = fragment.substr(0, m_minsearchlength);
	vector<Loc2pos> vec_of_poss = m_trieofDNA.find(tosearch, exactMatchOnly);
	
	for (int i = 0; i < vec_of_poss.size(); i++)
	{
		string segmentOfGenome;
		int locOfGenomeinArray = vec_of_poss[i].m_loc;
		int posOfFragmentinGenome = vec_of_poss[i].m_pos;
		int sizeOfRestOfGenome = m_genomes[locOfGenomeinArray].length() - posOfFragmentinGenome;
		
		for (int addon = 0; addon < sizeOfRestOfGenome; addon++)
		{
			if (m_genomes[locOfGenomeinArray].extract(posOfFragmentinGenome, minimumLength + addon, segmentOfGenome))
			{
				if (determineifExactMatchorSNiP(fragment.substr(0, segmentOfGenome.length()),segmentOfGenome) == 1)
				{
					DNAMatch newlymatched;
					newlymatched.genomeName = m_genomes[locOfGenomeinArray].name();
					newlymatched.length = segmentOfGenome.size();
					newlymatched.position = posOfFragmentinGenome;
					matches.push_back(newlymatched);
				}
				else if (determineifExactMatchorSNiP(fragment.substr(0, segmentOfGenome.length()), segmentOfGenome) == 2
					&& !exactMatchOnly)
				{
					DNAMatch newlymatched;
					newlymatched.genomeName = m_genomes[locOfGenomeinArray].name();
					newlymatched.length = segmentOfGenome.size();
					newlymatched.position = posOfFragmentinGenome;
					matches.push_back(newlymatched);
				}

			}
		}
	}
	
	if (matches.empty()) return false;
	return true;
}


bool GenomeMatcherImpl::findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const
{
	return false;  // This compiles, but may not be correct
}

//******************** GenomeMatcher functions ********************************

// These functions simply delegate to GenomeMatcherImpl's functions.
// You probably don't want to change any of this code.

GenomeMatcher::GenomeMatcher(int minSearchLength)
{
	m_impl = new GenomeMatcherImpl(minSearchLength);
}

GenomeMatcher::~GenomeMatcher()
{
	delete m_impl;
}

void GenomeMatcher::addGenome(const Genome& genome)
{
	m_impl->addGenome(genome);
}

int GenomeMatcher::minimumSearchLength() const
{
	return m_impl->minimumSearchLength();
}

bool GenomeMatcher::findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const
{
	return m_impl->findGenomesWithThisDNA(fragment, minimumLength, exactMatchOnly, matches);
}

bool GenomeMatcher::findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const
{
	return m_impl->findRelatedGenomes(query, fragmentMatchLength, exactMatchOnly, matchPercentThreshold, results);
}