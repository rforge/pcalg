/**
 * $Id: $
 */

#ifndef CONSTRAINT_HPP_
#define CONSTRAINT_HPP_

#include "constraint.hpp"
#include "gies_debug.hpp"

#include <boost/math/special_functions/log1p.hpp>
#include <boost/math/distributions/normal.hpp>

double IndepTestRFunction::test(uint u, uint v, std::vector<uint> S) const
{
	// Adapt indices to R convention
	std::vector<uint> shiftS;
	shiftS.reserve(S.size());
	std::vector<uint>::iterator vi;
	for (vi = S.begin(); vi != S.end(); ++vi)
		shiftS.push_back(*vi + 1);

	// Call R function to perform independence test
	return Rcpp::as<double>(_testFunction(u + 1, v + 1, shiftS, *_suffStat));
}

double IndepTestGauss::test(uint u, uint v, std::vector<uint> S) const
{
	// Calculate (absolute value of) z statistic
	#define CUT_THR 0.9999999
	double r, absz;
	//dout.level(3) << " Performing independence test for conditioning set of size " << S.size() << std::endl;
	if (S.empty())
		r = _correlation(u, v);
	else if (S.size() == 1) {
		uint w = S[0];
		r = (_correlation(u, v) - _correlation(u, w) * _correlation(v, w))/sqrt((1 - _correlation(v, w)*_correlation(v, w)) * (1 - _correlation(u, w)*_correlation(u, w)));
	}
	else {
		arma::mat PM;
		arma::uvec ind(S.size() + 2);
		ind(0) = u;
		ind(1) = v;
		for (uint i = 0; i < S.size(); ++i) ind(i + 2) = S[i];
		pinv(PM, _correlation.submat(ind, ind));
		// TODO include error handling
		r = - PM(0, 1)/sqrt(PM(0, 0) * PM(1, 1));
	}
	// Absolute value of r, respect cut threshold
	r = std::min(CUT_THR, std::abs(r));

	// Absolute value of z statistic
	// Note: log1p for more numerical stability, see "Aaux.R"; log1p is also available in
	// header <cmath>, but probably only on quite up to date headers (C++11)?
	absz = sqrt(_sampleSize - S.size() - 3) * 0.5 * boost::math::log1p(2*r/(1 - r));

	// Calculate p-value to z statistic (based on standard normal distribution)
	boost::math::normal distN;
	return (2*boost::math::cdf(boost::math::complement(distN, absz)));
}

bool Skeleton::hasEdge(const uint a, const uint b) const
{
	bool result;
	boost::tie(boost::tuples::ignore, result) = boost::edge(a, b, _graph);
	return result;
}


std::set<uint> Skeleton::getNeighbors(const uint vertex) const
{
	std::set<uint> result;
	UndirOutEdgeIter outIter, outLast;

	for (boost::tie(outIter, outLast) = boost::out_edges(vertex, _graph); outIter != outLast; outIter++)
			result.insert(boost::target(*outIter, _graph));

	return result;
}

Rcpp::LogicalMatrix Skeleton::getAdjacencyMatrix()
{
	Rcpp::LogicalMatrix result(getVertexCount(), getVertexCount());
	UndirEdgeIter ei, eiLast;
	for (boost::tie(ei, eiLast) = boost::edges(_graph); ei != eiLast; ei++) {
		dout.level(3) << "  Edge {" << boost::source(*ei, _graph) <<
				", " << boost::target(*ei, _graph) << "}\n";
		result(boost::source(*ei, _graph), boost::target(*ei, _graph)) = true;
		result(boost::target(*ei, _graph), boost::source(*ei, _graph)) = true;
	}

	return result;
}

void Skeleton::fitCondInd(const double alpha, int maxCondSize) {
	if (maxCondSize < 0)
		maxCondSize = getVertexCount();

	dout.level(2) << "Significance level " << alpha << std::endl;
	uint condSize, u, v, vtmp;
	int i;
	double pval;
	bool found = true;
	bool done = false;
	UndirOutEdgeIter outIter, outLast;
	std::set<uint>::iterator vi;
	std::vector<std::vector<uint>::iterator> si;
	std::vector<uint> condSet, neighbors;
	std::vector<uint>::iterator vnext;

	Skeleton oldSkel(getVertexCount());

	for (condSize = 1; found && condSize <= maxCondSize; ++condSize) {
		dout.level(2) << "Order = " << condSize << "; remaining edges: " << getEdgeCount() << std::endl;
		oldSkel._graph = _graph;
		condSet.resize(condSize);
		si.resize(condSize);
		found = false;
		for (u = 0; u < getVertexCount(); ++u)
			if (oldSkel.getDegree(u) > condSize) {
				// Found conditioning sets of size "condSize"
				found = true;
				dout.level(3) << "  Considering vertex " << u << " with degree " <<
						oldSkel.getDegree(u) << std::endl;

				// Initialize: first neighbor -> v, remaining neighbors -> vector of neighbors
				boost::tie(outIter, outLast) = boost::out_edges(u, oldSkel._graph);
				v = boost::target(*outIter, oldSkel._graph);
				neighbors.clear();
				neighbors.reserve(oldSkel.getDegree(u) - 1);
				for (outIter++; outIter != outLast; outIter++)
					neighbors.push_back(boost::target(*outIter, oldSkel._graph));
				vnext = neighbors.begin();

				// dout.level(3) << "  Neighbors: " << neighbors << std::endl;

				// Iterate over all neighbors of u
				done = false;
				do {
					// Initialize first conditioning set
					for (i = 0; i < condSize; ++i)
						si[i] = neighbors.begin() + i;

					// Iterate over conditioning sets
					do {
						for (i = 0; i < condSize; ++i)
							condSet[i] = *(si[i]);

						// Test of u and v are conditionally independent given condSet
						pval = _indepTest->test(u, v, condSet);
						dout.level(2) << "  x = " << u << ", y = " << v << ", S = " <<
								condSet << " : pval = " << pval << std::endl;
						if (pval >= alpha) {
							removeEdge(u, v);
							break; // Leave do-while-loop
						}

						// Proceed to next conditioning set
						for (i = condSize - 1;
								i >= 0 && si[i] == neighbors.begin() + (neighbors.size() - condSize + i);
								--i);
						if (i >= 0) {
							si[i]++;
							for (i++; i < condSize; ++i)
								si[i] = si[i - 1] + 1;
						}
					} while(i >= 0);

					// Proceed to next actual neighbor (note: "neighbors" is list of neighbors
					// from previous step)
					if (vnext != neighbors.end())
						do {
							vtmp = v;
							v = *vnext;
							*vnext = vtmp;
							vnext++;
						} while (!hasEdge(u, v) && vnext != neighbors.end());
					else
						done = true;
					// dout.level(3) << "  Neighbors: " << neighbors << std::endl;
				} while (!done);
			}
	}
}

#endif /* CONSTRAINT_HPP_ */
