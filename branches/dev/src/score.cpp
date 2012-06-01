/*
 * @author Alain Hauser
 * $Id: score.cpp 30 2012-02-26 15:04:56Z alhauser $
 */

#include "score.hpp"
#include "greedy.hpp"

#include <boost/tuple/tuple.hpp>
#include <limits>


double BICScoreRFunction::calcPartial(const uint vertex, std::set<uint> parents) const
{
	// Adapt indices to R convention...
	std::vector<uint> shiftParents;
	shiftParents.reserve(parents.size());
	std::set<uint>::iterator si;
	for (si = parents.begin(); si != parents.end(); ++si)
		shiftParents.push_back(*si + 1);

	// Call R function for partial score
	return Rcpp::as<double>(_partialScoreFunction(vertex + 1, shiftParents));
}

std::vector< std::vector<double> > BICScoreRFunction::calculateMLE(const EssentialGraph& dag) const
{
	std::vector< std::vector<double> > result(_vertexCount);

	uint v;
	std::set<uint> parents;
	std::set<uint>::iterator si;
	std::vector<uint> shiftParents;
	for (v = 0; v < _vertexCount; ++v) {
		// Get parents of vertex v and adapt their indices to the R convention
		parents = dag.getParents(v);
		shiftParents.clear();
		shiftParents.reserve(parents.size());
		for (si = parents.begin(); si != parents.end(); ++si)
			shiftParents.push_back(*si + 1);

		// Calculate parameters of vertex v using the R function
		result[v] = Rcpp::as< std::vector<double> >(_mleFunction(v + 1, shiftParents));
	}

	return result;
}

