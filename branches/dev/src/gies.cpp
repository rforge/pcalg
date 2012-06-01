/**
 * Main file of the Greedy Interventional Equivalence Search library for R
 *
 * @author Alain Hauser
 * $Id: gies.cpp 35 2012-03-11 13:34:35Z alhauser $
 */

#include "Rcpp.h"
#include <vector>
#include <string>
#include <algorithm>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

// For Mac g++ compiler... that bastard doesn't know uint...
typedef unsigned int uint;

#include "score.hpp"
#include "greedy.hpp"


using namespace boost::lambda;

/**
 * Extracts intervention targets from a SEXP to a TargetFamily object.
 */
TargetFamily castTargets(SEXP argTargets)
{
	int i;
	Rcpp::List listIventTargets(argTargets);
	std::vector<uint> vecTarget;
	std::vector<uint>::iterator vi;
	TargetFamily result;
	result.resize(listIventTargets.size());
	for (i = 0; i < listIventTargets.size(); i++) {
		vecTarget = listIventTargets[i];
		// Adapt indices to C++ convention...
		for (vi = vecTarget.begin(); vi != vecTarget.end(); ++vi)
			result[i].insert(*vi - 1);
	}
	return result;
}

/**
 * Reads in a graph from a list of in-edges passed as a SEXP to
 * an EssentialGraph object
 */
EssentialGraph castGraph(SEXP argInEdges)
{
	int i;
	Rcpp::List listInEdges(argInEdges);
	std::vector<uint> vecParents;
	std::vector<uint>::iterator vi;
	EssentialGraph result(listInEdges.size());
	for (i = 0; i < listInEdges.size(); ++i) {
		vecParents = listInEdges[i];
		// Adapt indices to C++ convention
		for (vi = vecParents.begin(); vi != vecParents.end(); ++vi)
			result.addEdge(*vi - 1, i);
	}

	return result;
}

/**
 * Wrap a graph structure to an R list of in-edges
 */
Rcpp::List wrapGraph(EssentialGraph graph)
{
	Rcpp::List result;
	Rcpp::IntegerVector vecEdges;
	std::set<uint> edges;
	std::set<uint>::iterator si;
	int i;

	for (i = 0; i < graph.getVertexCount(); ++i) {
		edges = graph.getInEdges(i);
		vecEdges = Rcpp::IntegerVector();
		for (si = edges.begin(); si != edges.end(); ++si)
			vecEdges.push_back(*si + 1);
		result.push_back(vecEdges);
	}

	return result;
}

/**
 * Interface to a variety of causal inference algorithms.
 *
 * @param 	argGraph			list of in-edges representing the current
 * 								essential graph
 * @param	argTargets			list of intervention targets, in the counting
 * 								convention of R: smallest vertex has index 1
 * @param 	argTargetIndex		list of intervention target indices, one per data
 * 								row. R counting convention
 * @param	argAlgorithm		string indicating the causal inference algorithm
 * 								to be used. Supported options: "GIES", "GDS", "DP"
 * @param	argScore			score objects to be used. Options: "scatter", "data"
 * 								for internal functions, or R function object
 * @param	argOptions			list of options specific for desired inference algorithm
 */
RcppExport SEXP causalInference(
		SEXP argGraph,
		SEXP argTargets,
		SEXP argAlgorithm,
		SEXP argScore,
		SEXP argOptions)
{
	// Initialize automatic exception handling; manual one does not work any more...
	initUncaughtExceptionHandler();

	// Cast graph
	EssentialGraph graph = castGraph(argGraph);

	// Cast list of targets
	TargetFamily targets = castTargets(argTargets);

	// Cast algorithm string and options
	std::string algName = Rcpp::as<std::string>(argAlgorithm);
	Rcpp::List options(argOptions);

	// TODO: cast score type, allow for C++ scoring objects
	// Up to now, only R functions are allowed for scoring...
	BICScore* score;
	score = new BICScoreRFunction(graph.getVertexCount(), &targets, Rcpp::Function(argScore), Rcpp::Function(argScore));

	graph.setScore(score);
	graph.setTargets(&targets);

	std::vector<int> steps;

	// Cast option for limits in vertex degree
	Rcpp::NumericVector maxDegree = options["maxdegree"];
	if (maxDegree.size() > 0) {
		if (maxDegree.size() == 1) {
			if (maxDegree[0] >= 1.) {
				uint uniformMaxDegree = static_cast<uint>(maxDegree[0]);
				graph.limitVertexDegree(uniformMaxDegree);
			}
			else {
				double maxRelativeDegree = maxDegree[0];
				graph.limitVertexDegree(maxRelativeDegree);
			}
		}
		else {
			std::vector<uint> maxDegrees = Rcpp::as< std::vector<uint> >(options["maxdegree"]);
			graph.limitVertexDegree(maxDegrees);
		}
	}
	
	// Cast option for vertices which are not allowed to have parents
	std::vector<uint> childrenOnly = Rcpp::as< std::vector<uint> >(options["childrenonly"]);
	std::for_each(childrenOnly.begin(), childrenOnly.end(), bind(&EssentialGraph::setChildrenOnly, &graph, _1, true));
	int stepLimit;

	// Perform inference algorithm:
	// GIES
	if (algName == "GIES") {
		// Enable caching, if requested
		if (options["caching"])
			graph.enableCaching();

		// Perform a greedy search, with or without turning phase
		// TODO: evtl. zusätzlichen Parameter einfügen, der wiederholtes Suchen
		// auch ohne Drehphase erlaubt...
		if (options["turning"]) {
			bool cont;
			do {
				cont = false;
				for (steps.push_back(0); graph.greedyForward(); steps.back()++);
				for (steps.push_back(0); graph.greedyBackward(); steps.back()++)
					cont = true;
				for (steps.push_back(0); graph.greedyTurn(); steps.back()++)
					cont = true;
			} while (cont);
		}
		else {
			for (steps.push_back(0); graph.greedyForward(); steps.back()++);
			for (steps.push_back(0); graph.greedyBackward(); steps.back()++);
		}
	}

	// Single phase or step of GIES
	else if (algName == "GIES-F" || algName == "GIES-B" || algName == "GIES-T") {
		// Limit to single step if requested
		stepLimit = options["maxsteps"];
		if (stepLimit == 0)
			stepLimit = graph.getVertexCount()*graph.getVertexCount();

		// Enable caching, if requested
		if (options["caching"])
			graph.enableCaching();

		steps.push_back(0);
		if (algName == "GIES-F")
			for (; steps.back() < stepLimit && graph.greedyForward(); steps.back()++);
		else if (algName == "GIES-B")
			for (; steps.back() < stepLimit && graph.greedyBackward(); steps.back()++);
		else if (algName == "GIES-T")
			for (; steps.back() < stepLimit && graph.greedyTurn(); steps.back()++);
	}

	// GDS
	else if (algName == "GDS") {
		// TODO: evtl. caching für GDS implementieren...
		// Perform a greedy search, with or without turning phase
		if (options["turning"]) {
			bool cont;
			do {
				cont = false;
				for (steps.push_back(0); graph.greedyDAGForward(); steps.back()++);
				for (steps.push_back(0); graph.greedyDAGBackward(); steps.back()++)
					cont = true;
				for (steps.push_back(0); graph.greedyDAGTurn(); steps.back()++)
					cont = true;
			} while (cont);
		}
		else {
			for (steps.push_back(0); graph.greedyDAGForward(); steps.back()++);
			for (steps.push_back(0); graph.greedyDAGBackward(); steps.back()++);
		}

		// Construct equivalence class
		graph.replaceUnprotected();
	}
	// DP
	else if (algName == "DP") {
		graph.dynamicProgrammingSearch();
		graph.replaceUnprotected();
	}

	// Return new list of in-edges and steps
	delete score;
	return Rcpp::List::create(Rcpp::Named("in.edges") = wrapGraph(graph),
			Rcpp::Named("steps") = steps);
}


//RcppExport SEXP optimalTarget(SEXP argAdjacency, SEXP argMaxSize)
//{
//	try {
//		// Cast adjacency matrix
//		Rcpp::IntegerMatrix matAdjacency(argAdjacency);
//		uint p = matAdjacency.nrow();
//		arma::imat adjacency(matAdjacency.begin(), matAdjacency.nrow(), matAdjacency.ncol(), false);
//
//		// Cast maximal size of intervention target
//		int maxSize = Rcpp::as<int>(argMaxSize);
//
//		// Create EssentialGraph instance
//		EssentialGraph graph(p);
//		graph.setAdjacencyMatrix(adjacency);
//
//		// Calculate optimal intervention target
//		std::set<uint> target = graph.getOptimalTarget(maxSize);
//
//		// Adapt numbering convention...
//		std::vector<uint> result(target.begin(), target.end());
//		std::for_each(result.begin(), result.end(), _1++);
//		return Rcpp::wrap(result);
//	} catch (std::exception& ex) {
//		forward_exception_to_r(ex);
//	}
//
//}
