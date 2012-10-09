/*
 * Classes for calculating the BIC score of essential graphs
 * given some data
 *
 * @author Alain Hauser
 * $Id$
 */

#ifndef SCORE_HPP_
#define SCORE_HPP_

#include "armaLapack.hpp"
#include <vector>
#include <set>
#include <boost/dynamic_bitset.hpp>
#include <Rcpp.h>

typedef unsigned int uint;

/**
 * Family of targets
 */
class TargetFamily : public std::vector<std::set<uint> >
{
public:
	/**
	 * Checks whether the family of targets protects a (hypothetical) arrow
	 * of a graph
	 */
	bool protects(const uint a, const uint b) const;
};
typedef TargetFamily::iterator TargetIterator;

/**
 * Extracts intervention targets from a SEXP to a TargetFamily object.
 */
TargetFamily castTargets(const SEXP argTargets);

/**
 * Casts a set of vertices (adapting the indexing convention) from a SEXP object.
 */
std::set<uint> castVertices(SEXP argVertices);

// Forward declaration for testing purpuses
class ScoreTest;

// Forward declaration
class EssentialGraph;

/**
 * Base scoring class
 */
class Score
{
	friend class ScoreTest;
protected:
	/**
	 * Number of random variables
	 */
	uint _vertexCount;

	/**
	 * Family of targets
	 */
	TargetFamily* _targets;
public:
	/**
	 * Constructors
	 */
	Score(uint vertexCount, TargetFamily* targets) :
		_vertexCount(vertexCount), _targets(targets) {};

	/**
	 * Virtual destructor
	 *
	 * Needed because of different storage of data.
	 */
	virtual ~Score() {}

	/**
	 * Virtual function yielding the total number of data points
	 */
	virtual uint getTotalDataCount() const
	{
		throw std::runtime_error("getTotalDataCount(): not implemented");
	}

	/**
	 * Virtual function yielding the number of data points
	 * available for estimating the conditional density of some given vertex
	 */
	virtual uint getDataCount(const uint vertex) const
	{
		throw std::runtime_error("getDataCount(uint): not implemented");
	}

	/**
	 * Virtual function to supply (preprocessed) data. Must be implemented by
	 * derived classes.
	 *
	 * @param 	data			preprocessed data as an R data structure.
	 * 							Exact type depends on type of score.
	 */
	virtual void setData(Rcpp::List& data) = 0;

	/**
	 * Calculates the local score of a vertex given its parents.
	 */
	virtual double local(const uint vertex, const std::set<uint>& parents) const = 0;

	/**
	 * Calculates the global score of a DAG.
	 */
	virtual double global(const EssentialGraph& dag) const;

	/**
	 * Calculates the local MLE for a vertex given its parents.
	 */
	virtual std::vector<double> localMLE(const uint vertex, const std::set<uint>& parents) const = 0;

	/**
	 * Calculates the MLE w.r.t. a DAG
	 *
	 * @param	dag		DAG with respect to which the MLE of the parameters are
	 * 					calculated
	 * @return			one vector of parameters per vertex of the DAG. The parameters
	 * 					describe the conditional density of the random variables given
	 * 					the values of their parents in the DAG. The concrete meaning of
	 * 					parameters is hence dependent on the chosen model class (i.e.,
	 * 					Gaussian, binary, ...)
	 */
	virtual std::vector< std::vector<double> > globalMLE(const EssentialGraph& dag) const = 0;
};

/**
 * Yields a pointer to a scoring object of appropriate type
 * @param	name		name of the scoring class; internal, short name used in R
 *  					library
 * @param	data		preprocessed data
 */
Score* createScore(std::string name, TargetFamily* targets, Rcpp::List data);

/**
 * Scoring class that acts as a wrapper to an external R function
 * doing the actual calculation
 */
class ScoreRFunction : public Score
{
protected:
	/**
	 * R function objects used to calculate: local score, global score, local MLE,
	 * global MLE
	 */
	Rcpp::Function _localScoreFunction;
	Rcpp::Function _globalScoreFunction;
	Rcpp::Function _localMLEFunction;
	Rcpp::Function _globalMLEFunction;
public:
	ScoreRFunction(uint vertexCount, TargetFamily* targets) :
			Score(vertexCount, targets) {};

	void setData(Rcpp::List& data);

	double local(const uint vertex, const std::set<uint>& parents) const;

	double global(const EssentialGraph& dag) const;

	std::vector<double> localMLE(const uint vertex, const std::set<uint>& parents) const;

	std::vector< std::vector<double> > globalMLE(const EssentialGraph& dag) const;
};

/**
 * Scoring class calculating a penalized l0-log-likelihood of Gaussian data,
 * based on precalculated scatter matrices.
 *
 * Special case: BIC score
 */
class ScoreGaussL0PenScatter : public Score
{
protected:
	/**
	 * Numbers of data points.
	 *
	 * For each vertex, the number of all data points coming from intervention NOT
	 * including this vertex are stored (n^{(i)}, 1 \leq i \leq p, in the usual
	 * notation).
	 */
	std::vector<int> _dataCount;
	uint _totalDataCount;

	/**
	 * Penalty constant
	 */
	double _lambda;

	/**
	 * Indicates whether an intercept should be calculated.
	 */
	bool _allowIntercept;

	/**
	 * Scatter matrices of interventional (or observational) data.
	 *
	 * The disjoint "cumulative" scatter matrices are stored, i.e. the sum of the scatter
	 * matrices of all interventions that do not contain a certain vertex.
	 * (In the usual notation, those are the matrices S^{(i)}, 1 \leq i \leq p)
	 */
	std::vector< arma::mat > _disjointScatterMatrices;

	/**
	 * Assignment of scatter matrices to vertices
	 */
	std::vector< arma::mat* > _scatterMatrices;

public:
	ScoreGaussL0PenScatter(uint vertexCount, TargetFamily* targets) :
		Score(vertexCount, targets),
		_dataCount(vertexCount),
		_scatterMatrices(vertexCount) {};

	uint getDataCount() const { return _totalDataCount;	}

	void setData(Rcpp::List& data);

	double local(const uint vertex, const std::set<uint>& parents) const;

	double global(const EssentialGraph& dag) const;

	std::vector<double> localMLE(const uint vertex, const std::set<uint>& parents) const;

	std::vector< std::vector<double> > globalMLE(const EssentialGraph& dag) const;
};

#endif /* SCORE_HPP_ */
