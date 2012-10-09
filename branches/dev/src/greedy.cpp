/**
 * greedy.cpp
 *
 * @author Alain Hauser
 * $Id$
 */

#include "greedy.hpp"
#include "gies_debug.hpp"

#include <algorithm>
#include <map>
#include <stack>
#include <boost/tuple/tuple.hpp>
#include <boost/unordered_map.hpp>

EssentialGraph::EssentialGraph(const uint vertexCount) :
	_graph(vertexCount),
	_maxVertexDegree(vertexCount, vertexCount),
	_childrenOnly(vertexCount)
{
	disableCaching();
}

void EssentialGraph::clear()
{
	boost::graph_traits<InternalEssentialGraph>::edge_iterator ei, ei_end, next;
	boost::tie(ei, ei_end) = boost::edges(_graph);
	for (next = ei; ei != ei_end; ei = next) {
		++next;
		boost::remove_edge(*ei, _graph);
	}
}

void EssentialGraph::addEdge(const uint a, const uint b, bool undirected)
{
	boost::add_edge(a, b, _graph);
	if (undirected)
		boost::add_edge(b, a, _graph);
}

void EssentialGraph::removeEdge(const uint a, const uint b, bool bothDirections)
{
	boost::remove_edge(a, b, _graph);
	if (bothDirections)
		boost::remove_edge(b, a, _graph);
}

bool EssentialGraph::existsPath(const uint a, const uint b, const std::set<uint>& C, const bool undirected)
{
	// Mark "forbidden" vertices as visited
	boost::dynamic_bitset<> visited(getVertexCount());
	std::set<uint>::iterator si;
	for (si = C.begin(); si != C.end(); ++si)
		visited.set(*si);

	// Trivial cases: if a or b are in C, return false
	if (visited.test(a) || visited.test(b))
		return false;

	// If (a, b) is an edge in the graph, remove it -- and add it again in the end
	bool restore = hasEdge(a, b);
	if (restore)
		removeEdge(a, b);

	// Check with depth-first search whether b is reachable from a without
	// using vertices in C
	std::stack<uint> nbhd;
	nbhd.push(a);
	visited.set(a);
	uint v;
	AdjacencyIter vi, vi_last;
	while (!nbhd.empty()) {
		v = nbhd.top();
		nbhd.pop();
		for (boost::tie(vi, vi_last) = boost::adjacent_vertices(v, _graph); vi != vi_last; ++vi)
			if (!undirected || hasEdge(*vi, v)) {
				if (*vi == b) {
					if (restore)
						addEdge(a, b);
					return true;
				}
				if (!visited.test(*vi)) {
					nbhd.push(*vi);
					visited.set(*vi);
				}
			}
	}

	if (restore)
		addEdge(a, b);
	return false;
}

bool EssentialGraph::existsPath(const uint a, const std::set<uint>& B, const std::set<uint>& C, const bool undirected) const
{
	// Exclude trivial cases
	if (B.empty() || std::includes(C.begin(), C.end(), B.begin(), B.end()) || C.find(a) != C.end())
		return false;

	// Mark "forbidden" vertices as visited
	boost::dynamic_bitset<> visited(getVertexCount());
	std::set<uint>::iterator si;
	for (si = C.begin(); si != C.end(); ++si)
		visited.set(*si);

	// Check with depth-first search whether any vertex in B is reachable from a without
	// using vertices in C
	std::stack<uint> nbhd;
	nbhd.push(a);
	visited.set(a);
	uint v;
	AdjacencyIter vi, vi_last;
	while (!nbhd.empty()) {
		v = nbhd.top();
		nbhd.pop();
		for (boost::tie(vi, vi_last) = boost::adjacent_vertices(v, _graph); vi != vi_last; ++vi)
			if (!undirected || hasEdge(*vi, v)) {
				if (B.find(*vi) != B.end())
					return true;
				if (!visited.test(*vi)) {
					nbhd.push(*vi);
					visited.set(*vi);
				}
			}
	}

	return false;
}

bool EssentialGraph::existsPath(const std::set<uint>& C, const uint a, const std::set<uint>& B)
{
	// Copy set of allowed vertices to bitset for faster lookup
	boost::dynamic_bitset<> allowed(getVertexCount());
	std::set<uint>::iterator si;
	for (si = C.begin(); si != C.end(); ++si)
		allowed.set(*si);

	// Exclude trivial cases
	std::set<uint> T = set_intersection(B, C);
	if (T.empty() || !allowed.test(a))
		return false;

	// Check with depth-first search whether any vertex in B is reachable from a
	// by only using vertices in C
	std::stack<uint> nbhd;
	boost::dynamic_bitset<> visited(getVertexCount());
	nbhd.push(a);
	visited.set(a);
	uint v;
	AdjacencyIter vi, vi_last;
	while (!nbhd.empty()) {
		v = nbhd.top();
		nbhd.pop();
		for (boost::tie(vi, vi_last) = boost::adjacent_vertices(v, _graph); vi != vi_last; ++vi)
			if (allowed.test(*vi)) {
				if (T.find(*vi) != T.end())
					return true;
				if (!visited.test(*vi)) {
					nbhd.push(*vi);
					visited.set(*vi);
				}
			}
	}

	return false;
}

std::vector<uint> EssentialGraph::greedyColoring(std::vector<uint> vertices)
{
	// Initialize coloring: vector of same length as list of vertices,
	// color of first vertex is 0
	std::vector<uint> coloring(vertices.size());
	boost::dynamic_bitset<> available;
	std::set<uint> adjacent;
	int i, j;

	for (i = 1; i < vertices.size(); ++i){
		// Assign vertex i the smallest color that has not yet been
		// used among its neighbors with smaller index
		adjacent = getAdjacent(vertices[i]);
		available.resize(adjacent.size());
		available.set();
		for (j = 0; j < i; ++j)
			if (isAdjacent(vertices[j], vertices[i]) && coloring[j] < adjacent.size())
				available.reset(coloring[j]);
		coloring[i] = (available.any() ? available.find_first() : adjacent.size());
	}

	return coloring;
}

bool EssentialGraph::hasEdge(const uint a, const uint b) const
{
	bool result;
	boost::tie(boost::tuples::ignore, result) = boost::edge(a, b, _graph);
	return result;
}

std::set<uint> EssentialGraph::getParents(const uint vertex) const
{
	std::set<uint> result;
	InEdgeIter ei, ei_end;
	Edge edge;

	for (boost::tie(ei, ei_end) = boost::in_edges(vertex, _graph); ei != ei_end; ei++) {
		edge = Edge(*ei, _graph);
		if (!hasEdge(edge.target, edge.source))
			result.insert(edge.source);
	}

	return result;
}

std::set<uint> EssentialGraph::getChildren(const uint vertex) const
{
	std::set<uint> result;
	OutEdgeIter ei, ei_end;
	Edge edge;

	for (boost::tie(ei, ei_end) = boost::out_edges(vertex, _graph); ei != ei_end; ei++) {
		edge = Edge(*ei, _graph);
		if (!hasEdge(edge.target, edge.source))
			result.insert(edge.target);
	}

	return result;
}

std::set<uint> EssentialGraph::getNeighbors(const uint vertex) const
{
	std::set<uint> result;
	InEdgeIter ei, ei_end;
	Edge edge;

	for (boost::tie(ei, ei_end) = boost::in_edges(vertex, _graph); ei != ei_end; ei++) {
		edge = Edge(*ei, _graph);
		if (hasEdge(edge.target, edge.source))
			result.insert(edge.source);
	}

	return result;
}

std::set<uint> EssentialGraph::getAdjacent(const uint vertex) const
{
	std::set<uint> result;
	InEdgeIter inIter, ei_end;
	OutEdgeIter outIter, outLast;

	for (boost::tie(inIter, ei_end) = boost::in_edges(vertex, _graph); inIter != ei_end; inIter++)
		result.insert(boost::source(*inIter, _graph));
	for (boost::tie(outIter, outLast) = boost::out_edges(vertex, _graph); outIter != outLast; outIter++)
		result.insert(boost::target(*outIter, _graph));

	return result;
}

std::set<uint> EssentialGraph::getInEdges(const uint vertex) const
{
	std::set<uint> result;
	InEdgeIter inIter, ei_end;

	for (boost::tie(inIter, ei_end) = boost::in_edges(vertex, _graph); inIter != ei_end; inIter++)
		result.insert(boost::source(*inIter, _graph));

	return result;
}

uint EssentialGraph::getDegree(const uint vertex) const
{
	return getAdjacent(vertex).size();
}

void EssentialGraph::limitVertexDegree(const std::vector<uint>& maxVertexDegree)
{
	if (maxVertexDegree.size() != getVertexCount())
		throw std::runtime_error("Number of vertex degrees must coincide with number of vertices");
	std::copy(maxVertexDegree.begin(), maxVertexDegree.end(), _maxVertexDegree.begin());
}

void EssentialGraph::limitVertexDegree(const uint maxVertexDegree)
{
	std::fill(_maxVertexDegree.begin(), _maxVertexDegree.end(), maxVertexDegree);
}

void EssentialGraph::limitVertexDegree(const double maxRelativeDegree)
{
	for (uint i = 0; i < getVertexCount(); ++i)
		_maxVertexDegree[i] = static_cast<uint>(maxRelativeDegree * (_score->getDataCount(i)));
}

boost::dynamic_bitset<> EssentialGraph::getAnteriorSet(const std::set<uint>& A)
{
	boost::dynamic_bitset<> result(getVertexCount());
    InEdgeIter ei, ei_end;
    Edge edge;
    std::set<uint>::iterator vi;
    std::stack<uint> nbhd;
    uint a;

    for (vi = A.begin(); vi != A.end(); ++vi) {
    	// Start a DFS at vertex *vi (one of the vertices in the start set)
    	nbhd.push(*vi);
    	result.set(*vi);
    	while (!nbhd.empty()) {
    		a = nbhd.top();
    		nbhd.pop();
    		// Move through the graph "backward", i.e. from edge targets to edge sources
    		for (boost::tie(ei, ei_end) = boost::in_edges(a, _graph); ei != ei_end; ++ei) {
    			edge = Edge(*ei, _graph);
    			// If newly detected source was not yet visited, add it to the neighborhood
    			// and to the result set (= anterior set)
    			if (!result.test(edge.source)) {
    				nbhd.push(edge.source);
    				result.set(edge.source);
    			}
    		}

    	}
    }

	return result;
}

boost::dynamic_bitset<> EssentialGraph::getPosteriorSet(const std::set<uint>& A)
{
	boost::dynamic_bitset<> result(getVertexCount());
    OutEdgeIter ei, ei_end;
    Edge edge;
    std::set<uint>::iterator vi;
    std::stack<uint> nbhd;
    uint a;

    for (vi = A.begin(); vi != A.end(); ++vi) {
    	// Start a DFS at vertex *vi (one of the vertices in the start set)
    	nbhd.push(*vi);
    	result.set(*vi);
    	while (!nbhd.empty()) {
    		a = nbhd.top();
    		nbhd.pop();
    		// Move through the graph "forward", i.e. from edge soruces to edge targets
    		for (boost::tie(ei, ei_end) = boost::out_edges(a, _graph); ei != ei_end; ++ei) {
    			edge = Edge(*ei, _graph);
    			// If newly detected target was not yet visited, add it to the neighborhood
    			// and to the result set (= posterior set)
    			if (!result.test(edge.target)) {
    				nbhd.push(edge.target);
    				result.set(edge.target);
    			}
    		}

    	}
    }

	return result;

}

//arma::umat EssentialGraph::getAdjacencyMatrix() const
//{
//	arma::umat adjacency(getVertexCount(), getVertexCount());
//	adjacency.zeros();
//	boost::graph_traits<InternalEssentialGraph>::edge_iterator ei, ei_end;
//	for (boost::tie(ei, ei_end) = boost::edges(_graph); ei != ei_end; ++ei)
//		adjacency(boost::source(*ei, _graph), boost::target(*ei, _graph)) = 1;
//
//	return adjacency;
//}

EssentialGraph EssentialGraph::getRepresentative() const
{
	EssentialGraph representative;
	representative._graph = _graph;

	// Orient all edges in the chain components according to a LexBFS-ordering
	boost::dynamic_bitset<> notVisited(getVertexCount());
	notVisited.set();
	std::set<uint> chainComp;
	std::set<uint>::iterator vi;
	uint v;
	while ((v = notVisited.find_first()) < getVertexCount()) {
		chainComp = representative.getChainComponent(v);
		representative.lexBFS(chainComp.begin(), chainComp.end(), true);
		for (vi = chainComp.begin(); vi != chainComp.end(); ++vi)
			notVisited.reset(*vi);
	}

	return representative;
}

void EssentialGraph::enableCaching()
{
	if (!_doCaching) {
		_doCaching = true;
		_actualPhase = CP_NONE;
		_scoreCache = std::vector<ArrowInsertion>(getVertexCount());
	}
}

void EssentialGraph::disableCaching()
{
	_doCaching = false;
	_actualPhase = CP_NONE;
	_scoreCache.clear();
}

std::set<uint> EssentialGraph::getChainComponent(const uint v) const
{
	std::vector<uint> nbhd(1, v);
	std::set<uint> chainComp;
	uint a;
	AdjacencyIter vi, vi_end;
	while (!nbhd.empty()) {
		a = nbhd.back();
		nbhd.pop_back();
		chainComp.insert(a);
		for (boost::tie(vi, vi_end) = boost::adjacent_vertices(a, _graph); vi != vi_end; vi++)
			if (hasEdge(*vi, a) && std::find(nbhd.begin(), nbhd.end(), *vi) == nbhd.end() && chainComp.find(*vi) == chainComp.end())
				nbhd.push_back(*vi);
	}
	return chainComp;
}

ArrowInsertion EssentialGraph::getOptimalArrowInsertion(const uint v)
{
	// For DEBUGGING purposes: print vertex being processed
	dout.level(2) << "Calculating optimal arrow insertion for vertex " << v << "\n";

	// Initialize optimal score gain
	ArrowInsertion result;
	result.score = 0.;

	// Respect maximum vertex degree and vertices that can only have children
	if (_childrenOnly.test(v) || (_maxVertexDegree[v] < getVertexCount() && getDegree(v) >= _maxVertexDegree[v]))
		return result;

	std::set<uint> C, C_par, C_sub, N;
	std::set<uint> neighbors, parents, adjacent;
	std::vector<std::set<uint> > maxCliques;
	std::set<uint>::iterator si;
	uint u;
	int i, j;
	double diffScore;
	CliqueStack cliqueStack;
	boost::dynamic_bitset<> posterior, forbidden;
	boost::unordered_map<std::set<uint>, double > localScore;
	boost::unordered_map<std::set<uint>, double >::iterator hmi;
	bool hasInserted;


	// Find maximal cliques in the neighborhood of v
	neighbors = getNeighbors(v);
	maxCliques = getMaxCliques(neighbors.begin(), neighbors.end());

	// Get parents of v (used for calculation of partial scores later on)
	parents = getParents(v);

	// Exclude forbidden sources:
	// - vertices reachable from children of v
	forbidden = getPosteriorSet(getChildren(v));
	// - vertices adjacent to v
	adjacent = getAdjacent(v);
	for (si = adjacent.begin(); si != adjacent.end(); ++si)
		forbidden.set(*si);
	// - v itself :-)
	forbidden.set(v);
	// - vertices which have reached the maximum degree
	for (u = 0; u < getVertexCount(); ++u)
		if (getDegree(u) >= _maxVertexDegree[u])
			forbidden.set(u);

	// Calculate vertices not reachable from v: for those, the "path condition"
	// for the clique C does not have to be checked later
	C.insert(v);
	posterior = getPosteriorSet(C);

	for (u = 0; u < getVertexCount(); ++u)
		if (!forbidden[u]) {
			// Calculate ne(v) \cap ad(u)
			N = set_intersection(neighbors, getAdjacent(u));

			// Add N as a set to check, and at the same time as a stop set.
			// Actually, N will be checked _only_ if it is a clique, i.e. subset
			// of a maximal clique
			cliqueStack.clear_all();
			cliqueStack.push_back(N);
			cliqueStack.stop_sets.insert(N);

			for (i = 0; i < maxCliques.size(); ++i) {
				// Only consider maximal cliques that contain N
				if (std::includes(maxCliques[i].begin(), maxCliques[i].end(), N.begin(), N.end())) {
					// Check all subsets of the actual maximal clique
					cliqueStack.append(maxCliques[i]);
					while(!cliqueStack.empty()) {
						C = cliqueStack.back();
						cliqueStack.pop_back();

						// Check whether there is a v-u-path that does not go through C
						if (!posterior.test(u) || !existsPath(v, u, C)) {
							// Calculate BIC score difference for current clique C
							// Note: if calculation is not possible (too low rank of
							// submatrices), local should return NaN; then the
							// test below fails
							// Use "localScore" as (additional) cache
							C_par = set_union(C, parents);
							hmi = localScore.find(C_par);
							if (hmi == localScore.end()) {
								dout.level(3) << "calculating partial score for vertex " << v << ", parents " << C_par << "...\n";
								diffScore = - _score->local(v, C_par);
								localScore[C_par] = diffScore;
							}
							else
								diffScore = hmi->second;
							dout.level(3) << "partial score for vertex " << v << ", parents " << C_par << ": " << -diffScore << "\n";
							C_par.insert(u);
							diffScore += _score->local(v, C_par);
							dout.level(3) << "partial score for vertex " << v << ", parents " << C_par << ": " << _score->local(v, C_par) << "\n";

							// If new score is better than previous optimum, and there is no
							// v-u-path that does not go through C, store (u, v, C) as new optimum
							if (diffScore > result.score) {
								result.source = u;
								result.clique = C;
								result.score = diffScore;
							}
						}

						// Add all subsets of C that differ from C in only one vertex to the stack
						for (si = C.begin(); si != C.end(); ++si) {
							C_sub = C;
							C_sub.erase(*si);
							cliqueStack.append(C_sub);
						}
					}

					// Mark the actual maximal set as checked, i.e. add it to the stop sets
					cliqueStack.stop_sets.insert(maxCliques[i]);
				}
			}
		}

	return result;
}

std::set<uint> EssentialGraph::_bitsToParents(const int vertex, const uint32_t bits)
{
	std::set<uint> parents;
	uint32_t pattern = 1;
	for (int i = 0; i < getVertexCount(); i++) {
		if (i != vertex) {
			if (bits & pattern)
				parents.insert(i);
			pattern *= 2;
		}
	}

	return parents;
}

std::vector<uint> EssentialGraph::_initialPEO(std::set<uint> chainComponent)
{

}

std::set<uint> EssentialGraph::_getOptimalUnrestrTarget()
{
	std::set<uint> target;

	// Get a greedy coloring of each chain component; return all vertices
	// with a color in the lower half
	boost::dynamic_bitset<> notVisited(getVertexCount());
	notVisited.set();
	std::set<uint> chainComp;
	std::vector<uint> ordering, coloring;
	std::vector<uint>::iterator vi;
	uint v, i, chromaticNumber;
	while ((v = notVisited.find_first()) < getVertexCount()) {
		// Get LexBFS-ordering of chain component of v
		chainComp = getChainComponent(v);
		ordering = lexBFS(chainComp.begin(), chainComp.end());

		// Get greedy coloring w.r.t. LexBFS-ordering
		coloring = greedyColoring(ordering);

		// Add vertices of lower half of coloring to the target
		chromaticNumber = *(std::max_element(coloring.begin(), coloring.end())) + 1;
		for (i = 0; i < ordering.size(); ++i) {
			if (coloring[i] < chromaticNumber/2)
				target.insert(ordering[i]);
			notVisited.reset(ordering[i]);
		}
	}

	return target;
}

uint EssentialGraph::_getOptimalSingleVertexTarget()
{
	uint u, v_opt, v, v_ind, i;
	uint eta, eta_opt, eta_min;
	std::set<uint> chainComp, chainCompVert, C, C_sub, neighbors;
	std::set<uint>::iterator si, sj;
	std::vector<uint> startOrder;
	std::vector<std::set<uint> > maxCliques;
	CliqueStack cliqueStack;
	EssentialGraph chainCompGraph, chainCompGraphCopy;
	TargetFamily targets;

	// Initialize optimal number of orientable edges
	eta_opt = 0;
	v_opt = getVertexCount();  // NO intervention helps anymore

	// Vertices that have not been checked yet
	boost::dynamic_bitset<> notChecked(getVertexCount());
	notChecked.set();

	// Search over all chain components
	while ((u = notChecked.find_first()) < getVertexCount()) {
		dout.level(3) << "  checking chain component of vertex " << u << "...\n";
		// Find chain component of v, and check whether it is non-trivial
		chainComp = getChainComponent(u);
		notChecked.reset(u);

		if (chainComp.size() > 1) {
			// Extract subgraph corresponding to chain component of v
			chainCompGraph = inducedSubgraph(chainComp.begin(), chainComp.end());

			v_ind = 0;
			for (si = chainComp.begin(); si != chainComp.end(); ++si, ++v_ind) {
				dout.level(3) << "  checking vertex " << *si << "...\n";
				notChecked.reset(*si);
				chainCompVert.clear();
				for (i = 0; i < chainComp.size(); ++i)
					chainCompVert.insert(i);

				// From now on, work in the induced subgraph

				// Set target family
				targets.clear();
				C.clear();
				targets.push_back(C);
				C.insert(v_ind);
				targets.push_back(C);
				chainCompGraph.setTargets(&targets);

				// Find maximal cliques in the neighborhood of v
				neighbors = chainCompGraph.getNeighbors(v_ind);
				maxCliques = chainCompGraph.getMaxCliques(neighbors.begin(), neighbors.end());

				// Initialize minimal number of orientable edges for given vertex
				eta_min = getVertexCount()*getVertexCount();

				// Check all cliques in the neighborhood of v...
				cliqueStack.clear_all();
				for (i = 0; i < maxCliques.size(); ++i) {
					cliqueStack.append(maxCliques[i]);
					while (!cliqueStack.empty()) {
						// ... and orient the edges of the chain component accordingly
						chainCompGraphCopy = chainCompGraph;
						C = cliqueStack.back();
						cliqueStack.pop_back();
						startOrder = std::vector<uint>(C.begin(), C.end());
						startOrder.push_back(v_ind);
						std::set_difference(chainCompVert.begin(), chainCompVert.end(), C.begin(), C.end(), std::inserter(startOrder, startOrder.end()));
						chainCompGraphCopy.lexBFS(startOrder.begin(), startOrder.end(), true);

						// Replace unprotected arrows (w.r.t. target v) in chain component
						chainCompGraphCopy.replaceUnprotected();
						eta = chainCompGraph.getEdgeCount() - chainCompGraphCopy.getEdgeCount();
						if (eta < eta_min)
							eta_min = eta;

						// Add all subsets of C that differ from C in only one vertex to the stack
						for (sj = C.begin(); sj != C.end(); ++sj) {
							C_sub = C;
							C_sub.erase(*sj);
							cliqueStack.append(C_sub);
						}
					} // WHILE !cliqueStack.empty()
				} // FOR i

				// Find maximal "eta_min"
				if (eta_min > eta_opt) {
					eta_opt = eta_min;
					v_opt = *si;
				}
			} // FOR si
		} // IF chainComp.size() > 1
	}

	return v_opt;
}

std::set<Edge, EdgeCmp> EssentialGraph::replaceUnprotected()
{
	// Map of arrow labels
	std::map<Edge, edge_flag, EdgeCmp> arrowFlags;
	// Set of undecidable arrows
	std::set<Edge, EdgeCmp> undecidableArrows;
	// Set of "affected" vertices, that is, vertices that are sources or targets
	// of an unprotected arrow
	std::set<Edge, EdgeCmp> result;

	Edge edge;

	// Find all arrows in the graph. Mark them as "protected", if they are
	// protected by an intervention target, and "undecidable" otherwise
	boost::graph_traits<InternalEssentialGraph>::edge_iterator edgeIter, edgeLast;
	for (boost::tie(edgeIter, edgeLast) = boost::edges(_graph); edgeIter != edgeLast; edgeIter++) {
		edge = Edge(*edgeIter, _graph);
		if (!hasEdge(edge.target, edge.source)) {
			if (_targets->protects(edge.source, edge.target)) arrowFlags[edge] = PROTECTED;
			else {
				undecidableArrows.insert(edge);
				arrowFlags[edge] = UNDECIDABLE;
			}
		}
	}

	// Check whether the arrows are part of a v-structure; if yes, mark them as "protected".
	std::map<Edge, edge_flag, EdgeCmp>::iterator arrIter1, arrIter2;
	uint v;
	for (v = 0; v < getVertexCount(); v++) {
		for (arrIter1 = arrowFlags.lower_bound(Edge(0, v)); arrIter1 != arrowFlags.upper_bound(Edge(getVertexCount(), v)); arrIter1++) {
			arrIter2 = arrIter1;
			for (arrIter2++; arrIter2 != arrowFlags.upper_bound(Edge(getVertexCount(), v)); arrIter2++)
				if (!isAdjacent((arrIter1->first).source, (arrIter2->first).source)) {
					arrIter1->second = PROTECTED;
					arrIter2->second = PROTECTED;
					undecidableArrows.erase(arrIter1->first);
					undecidableArrows.erase(arrIter2->first);
				}
		}
	}

	// Successively check all undecidable arrows, until no one remains
	std::set<Edge, EdgeCmp>::iterator undIter;
	edge_flag flag;
	while (!undecidableArrows.empty()) {
		// Find unprotected and protected arrows
		for (undIter = undecidableArrows.begin(); undIter != undecidableArrows.end(); undIter++) {
			edge = *undIter;
			flag = NOT_PROTECTED;

			// Check whether the arrow is in configuration (a)
			for (arrIter1 = arrowFlags.lower_bound(Edge(0, edge.source));
					flag != PROTECTED && arrIter1 != arrowFlags.upper_bound(Edge(getVertexCount(), edge.source));
					arrIter1++)
				if (!isAdjacent((arrIter1->first).source, edge.target))
					flag = (arrIter1->second == PROTECTED ? PROTECTED : UNDECIDABLE);

			// Check whether the arrow is in configuration (c)
			for (arrIter1 = arrowFlags.lower_bound(Edge(0, edge.target));
					flag != PROTECTED && arrIter1 != arrowFlags.upper_bound(Edge(getVertexCount(), edge.target));
					arrIter1++)
				if (isParent(edge.source, (arrIter1->first).source))
					flag = (arrIter1->second == PROTECTED && arrowFlags[Edge(edge.source, (arrIter1->first).source)] == PROTECTED ? PROTECTED : UNDECIDABLE);

			// Check whether the arrow is in configuration (d)
			for (arrIter1 = arrowFlags.lower_bound(Edge(0, edge.target));
					flag != PROTECTED && arrIter1 != arrowFlags.upper_bound(Edge(getVertexCount(), edge.target));
					arrIter1++) {
				arrIter2 = arrIter1;
				for (arrIter2++;
						flag != PROTECTED && arrIter2 != arrowFlags.upper_bound(Edge(getVertexCount(), edge.target));
						arrIter2++)
					if (isNeighbor(edge.source, (arrIter1->first).source) && isNeighbor(edge.source, (arrIter2->first).source) && !isAdjacent((arrIter1->first).source, (arrIter2->first).source))
						flag = (arrIter1->second == PROTECTED && arrIter2->second == PROTECTED ? PROTECTED : UNDECIDABLE);
			}

			// Store flag
			arrowFlags[edge] = flag;
		}

		// Replace unprotected arrows by lines; store affected edges in result set
		for (arrIter1 = arrowFlags.begin(); arrIter1 != arrowFlags.end(); ) {
			arrIter2 = arrIter1;
			arrIter1++;
			if (arrIter2->second != UNDECIDABLE)
				undecidableArrows.erase(arrIter2->first);
			if (arrIter2->second == NOT_PROTECTED) {
				addEdge((arrIter2->first).target, (arrIter2->first).source);
				result.insert(arrIter2->first);
				arrowFlags.erase(arrIter2);
			}
		}
	}

	return result;
}

void EssentialGraph::insert(const uint u, const uint v, const std::set<uint> C)
{
	// Temporary variables for caching
	std::set<Edge, EdgeCmp> directed, undirected, diffSet;
	std::set<Edge, EdgeCmp>::iterator ei;
	uint a;
	std::set<uint> recalc, recalcAnt;
	std::set<uint>::iterator si;
	boost::dynamic_bitset<> refreshCache(getVertexCount());
	EssentialGraph oldGraph;
	if (_doCaching)
		oldGraph = *this;

	// Get a LexBFS-ordering on the chain component of v, in which all edges of C
	// point toward v, and all other edges point away from v, and orient the edges
	// of the chain component accordingly
	std::set<uint> chainComp = getChainComponent(v);
	std::vector<uint> startOrder(C.begin(), C.end());
	startOrder.push_back(v);
	chainComp.erase(v);
	std::set_difference(chainComp.begin(), chainComp.end(), C.begin(), C.end(), std::inserter(startOrder, startOrder.end()));
	lexBFS(startOrder.begin(), startOrder.end(), true, &directed);

	// Add new arrow
	addEdge(u, v);

	// Successively replace unprotected arrows by lines
	undirected = replaceUnprotected();

	// If caching is enabled, recalculate the optimal arrow insertions where
	// necessary
	if (_doCaching) {
		// Genereate set of vertices whose anterior set is the set of vertices
		// whose cache has to be refreshed:
		// u, if there was no path from u to v before
		// TODO check conditions!!!
		// if (!oldGraph.existsPath(u, v))
			recalcAnt.insert(u);
		recalc.insert(u);
		// v, if the arrow was undirected and there was no path from v to u before
		// TODO check conditions!!
		// if (hasEdge(v, u) && !oldGraph.existsPath(v, u))
		if (hasEdge(v, u))
			recalcAnt.insert(v);
		recalc.insert(v);
		// the target of any newly directed edge
		diffSet = set_difference(directed, undirected);
		for (ei = diffSet.begin(); ei != diffSet.end(); ++ei) {
			recalcAnt.insert(ei->target);
			recalc.insert(ei->source);
		}
		// the source of any newly undirected edge
		diffSet = set_difference(undirected, directed);
		for (ei = diffSet.begin(); ei != diffSet.end(); ++ei) {
			recalcAnt.insert(ei->source);
			recalc.insert(ei->target);
		}

		// Calculate anterior set of that candidate set, and add vertices that
		// have to be recalculated without the complete anterior set
		refreshCache = getAnteriorSet(recalcAnt);
		for (si = recalc.begin(); si != recalc.end(); ++si)
			refreshCache.set(*si);

		// If v or u have reached the maximum degree, recalculate the optimal
		// arrow insertion for all vertices for which an insertion with new
		// parent u or v is proposed by the cache
		if (getDegree(u) >= _maxVertexDegree[u])
			for (a = 0; a < getVertexCount(); ++a)
				if (_scoreCache[a].source == u)
					refreshCache.set(a);
		if (getDegree(v) >= _maxVertexDegree[v])
			for (a = 0; a < getVertexCount(); ++a)
				if (_scoreCache[a].source == v)
					refreshCache.set(a);

		// Refresh cache: recalculate arrow insertions
		for (a = refreshCache.find_first(); a < getVertexCount(); a = refreshCache.find_next(a))
			_scoreCache[a] = getOptimalArrowInsertion(a);
	}
}

void EssentialGraph::remove(const uint u, const uint v, const std::set<uint> C)
{
	// Get a LexBFS-ordering on the chain component of v, in which all edges of C
	// (or C \cup \{u\}, if u lies in the chain component)
	// point toward v, and all other edges point away from v, and orient the edges
	// of the chain component accordingly
	std::set<uint> chainComp = getChainComponent(v);
	std::vector<uint> startOrder(C.begin(), C.end());
	if (chainComp.find(u) != chainComp.end())
		startOrder.push_back(u);
	startOrder.push_back(v);
	chainComp.erase(v);
	std::set_difference(chainComp.begin(), chainComp.end(), C.begin(), C.end(), std::inserter(startOrder, startOrder.end()));
	lexBFS(startOrder.begin(), startOrder.end(), true);

	// Remove the edge between u and v
	removeEdge(u, v, true);

	// Successively replace unprotected arrows by lines
	replaceUnprotected();
}

void EssentialGraph::turn(const uint u, const uint v, const std::set<uint> C)
{
	std::set<uint> chainComp;
	std::vector<uint> startOrder;

	// If u and v lie in different chain components, first order the edges in the
	// chain component of u
	if (!hasEdge(u, v)) {
		chainComp = getChainComponent(u);
		chainComp.erase(u);
		startOrder.push_back(u);
		startOrder.insert(startOrder.end(), chainComp.begin(), chainComp.end());
		lexBFS(startOrder.begin(), startOrder.end(), true);
		startOrder.clear();
	}

	// Get a LexBFS-ordering on the chain component of v in which all edges of C
	// point towards v, and all other edges point away from v
	chainComp = getChainComponent(v);
	startOrder.insert(startOrder.end(), C.begin(), C.end());
	startOrder.push_back(v);
	chainComp.erase(v);
	if (hasEdge(u, v)) {
		startOrder.push_back(u);
		chainComp.erase(u);
	}
	std::set_difference(chainComp.begin(), chainComp.end(), C.begin(), C.end(), std::inserter(startOrder, startOrder.end()));
	lexBFS(startOrder.begin(), startOrder.end(), true);

	// Turn the edge (v, u)
	removeEdge(v, u);
	addEdge(u, v);

	// Successively replace unprotected arrows by lines
	replaceUnprotected();
}

bool EssentialGraph::greedyForward()
{
	uint v, v_opt;
	std::vector<ArrowInsertion>::iterator si;
	ArrowInsertionCmp comp;
	ArrowInsertion insertion, optInsertion;

	// For DEBUGGING purposes: print phase
	dout.level(3) << "== starting forward phase...\n";

	// Initialize optimal score gain
	optInsertion.score = 0.;

	// If caching is disabled calculate the score differences for all possible edges
	if (!_doCaching) {
		for (v = 0; v < getVertexCount(); v++) {
			// Calculate optimal arrow insertion for given target vertex v
			insertion = getOptimalArrowInsertion(v);

			// Look for optimal score
			if (insertion.score > optInsertion.score) {
				optInsertion = insertion;
				v_opt = v;
			}
		}
	}

	// If caching is enabled, search the cache for the best score
	dout.level(3) << "vertex count: " << getVertexCount() << "\n";
	if (_doCaching) {
		// If score has to be initialized (current phase is not forward), do it
		if (_actualPhase != CP_FORWARD)
			for (v = 0; v < getVertexCount(); v++)
				_scoreCache[v] = getOptimalArrowInsertion(v);

		// Find optimal arrow insertion from cache
		si = std::max_element(_scoreCache.begin(), _scoreCache.end(), comp);
		v_opt = std::distance(_scoreCache.begin(), si);
		optInsertion = *si;
		_actualPhase = CP_FORWARD;
	}

	// If the score can be augmented, do it
	if (optInsertion.score > 0.) {
		// For DEBUGGING purposes: print inserted arrow
		dout.level(3) << "  inserting edge (" << optInsertion.source << ", " << v_opt << ") with C = "
				<< optInsertion.clique << ", S = " << optInsertion.score << "\n";
		insert(optInsertion.source, v_opt, optInsertion.clique);
//		#if DEBUG_OUTPUT_LEVEL >= 2
//			getAdjacencyMatrix().print("A = ");
//		#endif
		return true;
	}
	else
		return false;
}

bool EssentialGraph::greedyBackward()
{
	std::set<uint> neighbors, parents, candidates;
	std::vector<std::set<uint> > maxCliques;
	std::set<uint> C, C_par, C_sub, C_opt, N;
	std::set<uint>::iterator iter, ui;
	int i, j;
	uint v, u_opt, v_opt;
	double diffScore, diffScore_opt;
	CliqueStack cliqueStack;

	// For DEBUGGING purposes: print phase
	dout.level(3) << "== starting backward phase...\n" ;

	// Initialize optimal score gain
	diffScore_opt = 0.;

	for (v = 0; v < getVertexCount(); v++) {
		// Get neighbors and parents of v (used for calculation of BIC score later on)
		neighbors = getNeighbors(v);
		parents = getParents(v);
		candidates = set_union(neighbors, parents);
		for (ui = candidates.begin(); ui != candidates.end(); ++ui) {
			N = set_intersection(neighbors, getAdjacent(*ui));

			// Find all maximal cliques on N
			maxCliques = getMaxCliques(N.begin(), N.end());
			cliqueStack.clear_all();

			// Calculate the score difference for all cliques in N
			for (i = 0; i < maxCliques.size(); ++i) {
				// Check all subsets of the actual maximal clique
				cliqueStack.append(maxCliques[i]);
				while(!cliqueStack.empty()) {
					C = cliqueStack.back();
					cliqueStack.pop_back();

					// Calculate BIC score difference for current clique C
					C_par = set_union(C, parents);
					C_par.insert(*ui);
					diffScore = - _score->local(v, C_par);
					C_par.erase(*ui);
					diffScore += _score->local(v, C_par);

					// If new score is better than previous optimum, store (u, v, C) as new optimum
					if (diffScore > diffScore_opt) {
						u_opt = *ui;
						v_opt = v;
						C_opt = C;
						diffScore_opt = diffScore;
					}

					// Add all subsets of C that differ from C in only one vertex to the stack
					for (iter = C.begin(); iter != C.end(); ++iter) {
						C_sub = C;
						C_sub.erase(*iter);
						cliqueStack.append(C_sub);
					}
				} // while !cliqueStack.empty()

				// Mark the actual maximal set as checked, i.e. add it to the stop sets
				cliqueStack.stop_sets.insert(maxCliques[i]);
			} // for i
		} // for ui
	} // for v

	// If caching is enabled, store current phase...
	// TODO: Change that to admit actual caching also in the backward phase!!!
	if (_doCaching)
		_actualPhase = CP_BACKWARD;

	// If the score can be augmented, do it
	if (diffScore_opt > 0.) {
		// For DEBUGGING purposes
		dout.level(1) << "  deleting edge (" << u_opt << ", " << v_opt << ") with C = "
				<< C_opt << ", S = " << diffScore_opt << "\n";
		remove(u_opt, v_opt, C_opt);
		//getAdjacencyMatrix().print("A = ");
		return true;
	}
	else
		return false;
}

bool EssentialGraph::greedyTurn()
{
	std::set<uint> children, neighbors, neighbors2, parents, C, C_par, C_opt, C_sub, CminN, N;
	std::vector<std::set<uint> > maxCliques;
	std::set<uint>::iterator iter, ui;
	uint v, u_opt, v_opt;
	int i;
	double diffScore, diffScore_opt;
	CliqueStack cliqueStack;

	// For DEBUGGING purposes: print phase
	dout.level(3) << "== starting turning phase...\n";

	// Initialize optimal score gain
	diffScore_opt = 0.;

	for (v = 0; v < getVertexCount(); ++v) 
		// Respect vertices that are not allowed to have parents, but only children
		if (!_childrenOnly.test(v))  {
			dout.level(3) << "  checking edges incident to v = " << v << "...\n";

			// Store parents and neighbors of v
			neighbors = getNeighbors(v);
			parents = getParents(v);

			// Search over all neighbors of v; turning non-essential arrows
			for (ui = neighbors.begin(); ui != neighbors.end(); ++ui) {
				N = set_intersection(neighbors, getNeighbors(*ui));

				// Find all maximal cliques in the neighborhood of v (without u)
				neighbors2 = neighbors;
				neighbors2.erase(*ui);
				maxCliques = getMaxCliques(neighbors2.begin(), neighbors2.end());
				cliqueStack.clear_all();

				// Calculate the score difference for all (admissible) cliques in the neighborhood of v
				for (i = 0; i < maxCliques.size(); ++i) {
					// Check all subsets of the actual maximal clique
					cliqueStack.append(maxCliques[i]);
					while(!cliqueStack.empty()) {
						C = cliqueStack.back();
						cliqueStack.pop_back();

						// Check if C \ N is not empty, and if C \cap N separates C \ N and N \ C in
						// the neighborhood of v
						CminN = set_difference(C, N);
						if (!(CminN.empty()) && !existsPath(set_difference(neighbors, set_intersection(C, N)), *(CminN.begin()), set_difference(N, C))) {
							// Calculate BIC score difference for current clique C
							C_par = set_union(C, parents);
							diffScore = - _score->local(v, C_par);
							C_par.insert(*ui);
							diffScore += _score->local(v, C_par);
							C_par = set_union(set_intersection(C, N), getParents(*ui));
							diffScore += _score->local(*ui, C_par);
							C_par.insert(v);
							diffScore -= _score->local(*ui, C_par);

							// If new score is better than previous optimum, store (u, v, C) as new optimum
							if (diffScore > diffScore_opt) {
								u_opt = *ui;
								v_opt = v;
								C_opt = C;
								diffScore_opt = diffScore;
							}
						}

						// Add all subsets of C that differ from C in only one vertex to the stack
						for (iter = C.begin(); iter != C.end(); ++iter) {
							C_sub = C;
							C_sub.erase(*iter);
							cliqueStack.append(C_sub);
						}
					} // while !cliqueStack.empty()

					// Mark the actual maximal set as checked, i.e. add it to the stop sets
					cliqueStack.stop_sets.insert(maxCliques[i]);
				} // for i
			} // for ui

			// Find all maximal cliques in the neighborhood of v
			maxCliques = getMaxCliques(neighbors.begin(), neighbors.end());

			// Search over all children of v; turning essential arrows
			children = getChildren(v);
			for (ui = children.begin(); ui != children.end(); ++ui) {
				dout.level(2) << "  trying to turn arrow (" << v << ", " << *ui << ")...\n";
				N = set_intersection(neighbors, getParents(*ui));

				// Add N as a set to check, and at the same time as a stop set.
				// Actually, N will be checked _only_ if it is a clique, i.e. subset
				// of a maximal clique
				cliqueStack.clear_all();
				cliqueStack.push_back(N);
				cliqueStack.stop_sets.insert(N);

				for (i = 0; i < maxCliques.size(); ++i) {
					// Only consider maximal cliques that contain N
					if (std::includes(maxCliques[i].begin(), maxCliques[i].end(), N.begin(), N.end())) {
						// Check all subsets of the actual maximal clique
						cliqueStack.append(maxCliques[i]);
						while(!cliqueStack.empty()) {
							C = cliqueStack.back();
							cliqueStack.pop_back();

							// Check if there is no v-u-path (except (v, u)) that does not visit C \cup ne(u)
							if (!existsPath(v, *ui, set_union(C, getNeighbors(*ui)))) {
								// Calculate BIC score difference for current clique C
								C_par = set_union(C, parents);
								diffScore = - _score->local(v, C_par);
								C_par.insert(*ui);
								diffScore += _score->local(v, C_par);
								C_par = getParents(*ui);
								diffScore -= _score->local(*ui, C_par);
								C_par.erase(v);
								diffScore += _score->local(*ui, C_par);
								dout.level(2) << "  score difference for (u, v, C) = (" << *ui << ", " << v << ", " << C << "): " << diffScore << "\n";


								// If new score is better than previous optimum, and there is no
								// v-u-path that does not go through C, store (u, v, C) as new optimum
								if (diffScore > diffScore_opt) {
									u_opt = *ui;
									v_opt = v;
									C_opt = C;
									diffScore_opt = diffScore;
								}
							}

							// Add all subsets of C that differ from C in only one vertex to the stack
							for (iter = C.begin(); iter != C.end(); ++iter) {
								C_sub = C;
								C_sub.erase(*iter);
								cliqueStack.append(C_sub);
							}
						}

						// Mark the actual maximal set as checked, i.e. add it to the stop sets
						cliqueStack.stop_sets.insert(maxCliques[i]);
					}
				}
			} // for ui
		}

	// If caching is enabled, store current phase...
	// TODO: Change that to admit actual caching also in the turning phase!!!
	if (_doCaching)
		_actualPhase = CP_TURNING;

	// Turn the highest-scoring edge
	if (diffScore_opt > 0.) {
		// For DEBUGGING purposes
		dout.level(1) << "  turning edge (" << v_opt << ", " << u_opt << ") with C = " << C_opt << ", S = " << diffScore_opt << "\n";
		turn(u_opt, v_opt, C_opt);
//		#if DEBUG_OUTPUT_LEVEL >= 2
//			getAdjacencyMatrix().print("A = ");
//		#endif
		return true;
	}
	else
		return false;

}

bool EssentialGraph::greedyDAGForward()
{
	uint u, v, u_opt, v_opt, p;
	double diffScore, diffScore_opt;
	std::set<uint> parents, C_new;

	// Initialize help variables
	diffScore_opt = 0.;
	p = getVertexCount();

	// Find edge that maximally increases BIC score when added
	for (v = 0; v < p; ++v) {
		parents = getParents(v);
		for (u = 0; u < p; ++u)
			if (u != v && !isAdjacent(u, v) && !existsPath(v, u)) {
				// Calculate BIC score difference for adding edge (u, v)
				C_new = parents;
				diffScore = - _score->local(v, C_new);
				C_new.insert(u);
				diffScore += _score->local(v, C_new);

				// If new score is better than previous optimum
				if (diffScore > diffScore_opt) {
					u_opt = u;
					v_opt = v;
					diffScore_opt = diffScore;
				}
			}
	}

	// Add this edge, if it incrases the BIC score
	if (diffScore_opt > 0.) {
		addEdge(u_opt, v_opt);
		return true;
	}
	else
		return false;
}

bool EssentialGraph::greedyDAGBackward()
{
	uint u, v, u_opt, v_opt, p;
	double diffScore, diffScore_opt;
	std::set<uint> parents, C_new;
	std::set<uint>::iterator ui;

	// Initialize help variables
	diffScore_opt = 0.;
	p = getVertexCount();

	// Find edge that maximally increases BIC score when removed
	for (v = 0; v < p; ++v) {
		parents = getParents(v);
		for (ui = parents.begin(); ui != parents.end(); ++ui) {
			// Calculate BIC score difference when removing edge (u, v)
			C_new = parents;
			diffScore = - _score->local(v, C_new);
			C_new.erase(*ui);
			diffScore += _score->local(v, C_new);

			// If new score is better than previous optimum, store (u, v, C) as new optimum
			if (diffScore > diffScore_opt) {
				u_opt = *ui;
				v_opt = v;
				diffScore_opt = diffScore;
			}
		}
	}

	// Remove this edge, if this incrases the BIC score
	if (diffScore_opt > 0.) {
		removeEdge(u_opt, v_opt);
		return true;
	}
	else
		return false;
}

bool EssentialGraph::greedyDAGTurn()
{
	uint v, u_opt, v_opt, p;
	double diffScore, diffScore_opt;
	std::set<uint> parents, C_new, D_new, emptyset;
	std::set<uint>::iterator ui;

	// Initialize help variables
	diffScore_opt = 0.;
	p = getVertexCount();

	// Find edge that maximally increases BIC score when turned, i.e. when its
	// orientation is changed
	for (v = 0; v < p; ++v) {
		parents = getParents(v);
		for (ui = parents.begin(); ui != parents.end(); ++ui) {
			if (!existsPath(*ui, v)) {
				C_new = parents;
				D_new = getParents(*ui);
				diffScore = - _score->local(v, C_new) - _score->local(*ui, D_new);
				C_new.erase(*ui);
				D_new.insert(v);
				diffScore += _score->local(v, C_new) + _score->local(*ui, D_new);

				// If new score is better than previous optimum, store (u, v, C) as new optimum
				if (diffScore > diffScore_opt) {
					u_opt = *ui;
					v_opt = v;
					diffScore_opt = diffScore;
				}
			}
		}
	}

	// Turn this edge, if this incrases the BIC score
	if (diffScore_opt > 0.) {
		removeEdge(u_opt, v_opt);
		addEdge(v_opt, u_opt);
		return true;
	}
	else
		return false;
}

void EssentialGraph::dynamicProgrammingSearch()
{
	// Check whether DAG is not too large (practically, the algorithm will
	// only work for substantially smaller graphs than p = 32 due to
	// memory and run time limits)
	if (getVertexCount() > 31)
		throw std::length_error("Vertex count must not exceed 31.");

	int i, j;
	uint32_t subset, subsubset, pattern;
	bool unset;

	// Table of best parents given a variable and a set of candidate parents
	std::vector< std::vector<uint32_t> > bestParents(getVertexCount(), std::vector<uint32_t>(1 << (getVertexCount() - 1)));
	// for (i = 0; i < bestParents.size(); i++)
	//	bestParents[i].resize(1 << (_vertexCount - 1), 0);
	// Table with corresponding best local scores
	std::vector< std::vector<double> > localScore(getVertexCount(), std::vector<double>(1 << (getVertexCount() - 1)));
	// for (i = 0; i < localScore.size(); i++)
	//	localScore[i].resize(1 << (_vertexCount - 1), 0.);

	// Table of best sinks for all subsets of variables
	std::vector< int > bestSink(1 << getVertexCount());
	// Table of corresponding optimal scores
	std::vector< double > totalScore(1 << getVertexCount(), 0.);

	// Forward phase of DP: fill in tables of best parents and corresponding local scores
	for (i = 0; i < getVertexCount(); i++) {
		//std::cout << "\n" << i << ":\t";

		for (subset = 0; subset < bestParents[i].size(); subset++) {
			localScore[i][subset] = _score->local(i, _bitsToParents(i, subset));
			//std::cout << "localscore: " << localScore[i][subset] << "\n";
			bestParents[i][subset] = subset;
			for (pattern = 1; pattern < bestParents[i].size(); pattern *= 2) {
				subsubset = subset & (~pattern);
				//std::cout << "subsubset: " << subsubset << "\n";
				if (localScore[i][subsubset] > localScore[i][subset]) {
					localScore[i][subset] = localScore[i][subsubset];
					bestParents[i][subset] = bestParents[i][subsubset];
				}
			}

			//std::cout << "(" << bestParents[i][subset] << ", " << localScore[i][subset] << ")\t";
		}
	}


	// Fill in tables of best sinks and corresponding total scores
	//std::cout << "\n\nBest sinks:\n";
	for (subset = 1; subset < totalScore.size(); subset++) {
		unset = true;
		for (i = 0; i < getVertexCount(); i++) {
			if ((1 << i) & subset) {
				// Calculate subsubset w.r.t. candidate set for parents, not w.r.t. full set of variables
				pattern = (1 << i) - 1;
				subsubset = (subset & pattern) | ((subset & (~pattern << 1)) >> 1);
				//std::cout << "subset: " << subset << ";\ti: " << i << ";\tpattern: " << pattern << ";\tsubsubset: " << subsubset << "\n";
				if (unset || totalScore[subset] < localScore[i][subsubset] + totalScore[subset & ~(1 << i)]) {
					totalScore[subset] = localScore[i][subsubset] + totalScore[subset & ~(1 << i)];
					bestSink[subset] = i;
					unset = false;
				}
			}
		}

		//std::cout << subset << ": " << bestSink[subset] << ", " << totalScore[subset] << "\n";
	}

	// Backward phase of DP: reconstruct optimal DAG
	clear();
	std::set<uint> parents;
	std::set<uint>::iterator pi;
	subset = (1 << getVertexCount()) - 1;
	while (subset) {
		i = bestSink[subset];

		// Calculate subsubset w.r.t. candidate set for parents, not w.r.t. full set of variables
		pattern = (1 << i) - 1;
		subsubset = (subset & pattern) | ((subset & (~pattern << 1)) >> 1);
		//std::cout << "subset: " << subset << ";\tbest sink: " << i << ";\tsubsubset: " << subsubset << "\n";

		parents = _bitsToParents(i, bestParents[i][subsubset]);
		for (pi = parents.begin(); pi != parents.end(); ++pi)
			addEdge(*pi, i);

		subset = subset & ~(1 << i);
	}
}


std::set<uint> EssentialGraph::getOptimalTarget(uint maxSize)
{
	std::set<uint> target;

	if (maxSize == 0)
		return target;
	else if (maxSize == 1) {
		uint v = _getOptimalSingleVertexTarget();
		if (v < getVertexCount())
			target.insert(v);
		return target;
	}
	else if (maxSize == getVertexCount())
		return _getOptimalUnrestrTarget();
	else
		throw std::runtime_error("Optimal targets with size other than 1 or p are not supported.");
}
