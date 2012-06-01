/*
 * Classes for calculating the BIC score of essential graphs
 * given some data
 *
 * @author Alain Hauser
 * $Id: score.hpp 35 2012-03-11 13:34:35Z alhauser $
 */

#ifndef SCORE_HPP_
#define SCORE_HPP_

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

// Forward declaration for testing purpuses
class BICScoreTest;

// Forward declaration
class EssentialGraph;

/**
 * Base BIC scoring class
 */
class BICScore
{
	friend class BICScoreTest;
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
	BICScore(uint vertexCount, TargetFamily* targets) :
		_vertexCount(vertexCount), _targets(targets) {};

	/**
	 * Virtual destructor
	 *
	 * Needed because of different storage of data.
	 */
	virtual ~BICScore() {}

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
	 * Virtual function to supply (interventional) data. Must be implemented by
	 * derived classes.
	 *
	 * @param	targetIndex		vector of intervention target indices. MUST BE
	 * 							IN ASCENDING ORDER!
	 * @param 	data			n x p matrix of data
	 */
	// NOTE: no virtual method any more since different types of score objects
	// may use different type of data matrices (e.g., double, integer, ...)
	// virtual void setData(const std::vector<uint>& targetIndex, arma::mat& data) = 0;

	/**
	 * Calculates the partial BIC score of a vertex given its parents.
	 */
	virtual double calcPartial(const uint vertex, std::set<uint> parents) const = 0;

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
	virtual std::vector< std::vector<double> > calculateMLE(const EssentialGraph& dag) const = 0;
};

// TODO: add classes for calculating the score in C++
// (to be copied from the old sources in the GIES library)

/**
 * BIC scoring class that acts as a wrapper to an external R function
 * doing the actual calculation
 */
class BICScoreRFunction : public BICScore
{
protected:
	/**
	 * Total number of data points
	 */
	// uint _totalDataCount;

	/**
	 * R function object used to calculate the partial score
	 */
	Rcpp::Function _partialScoreFunction;


	/**
	 * R function object used to calculate the MLE
	 */
	Rcpp::Function _mleFunction;
public:
	BICScoreRFunction(uint vertexCount,
		TargetFamily* targets,
		const Rcpp::Function partialScoreFunction,
		const Rcpp::Function mleFunction) :
		BICScore(vertexCount, targets),
		_partialScoreFunction(partialScoreFunction),
		_mleFunction(mleFunction) {};

	// uint getTotalDataCount() const { return _totalDataCount; }

	double calcPartial(const uint vertex, std::set<uint> parents) const;

	std::vector< std::vector<double> > calculateMLE(const EssentialGraph& dag) const;
};


#endif /* SCORE_HPP_ */
