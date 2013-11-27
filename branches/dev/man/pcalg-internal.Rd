\name{pcalg-internal}
%
%\alias{checkEdges} not exported
%\alias{lm.cov}
\alias{adj.check}
\alias{ancTS}
\alias{caus.inf}
\alias{check.new.coll}
\alias{checkTriple}
\alias{ci.test}
\alias{dep.triple}
\alias{dreach}%-> ./dreach.old
\alias{faith.check}
\alias{find.min.discr.path}
\alias{find.sink}
\alias{find.unsh.triple}
\alias{flipEdges}
\alias{legal.psep}
\alias{minUncovCircPath}
\alias{minUncovPdPath}
\alias{my.SpecialDag}
\alias{pc.cons.intern}%-. ./pc.cons.intern.old
\alias{pcAlgo.Perfect}
\alias{pdsep}%-> ./pdsep.old;  mentioned on ./fci.Rd
\alias{possibleDe}
\alias{qreach}%-> ./qreach.old; mentioned on ./fci.Rd
\alias{reach}
\alias{rfci.vstructures}
\alias{triple2numb}
\alias{udag2apag}%-> ./udag2apag.old
\alias{udag2pag} %-> ./udag2pag.old
\alias{udag2pdagRelaxed}%-> ./udag2pdagRelaxed.old
\alias{udag2pdagSpecial}%-> ./udag2pdagSpecial.old
\alias{updateList}
\alias{visibleEdge}
\alias{wgtMatrix}
%--------------------------------------------------------------------------
% FIXME: *do* document [says Martin]  (((or confirm current status))):
\alias{gAlgo-class}
\alias{trueCov}% used in examples / tests!
\alias{allDags}%<- mentioned in other help pages !!
\alias{amat2dag}% in ../tests/test_amat2dag.R
%% the following seem important, self contained, have "header help"
\alias{skeleton.dag2pag}
%
\title{Internal Pcalg Functions}
\description{
  Internal functions of package \pkg{pcalg}
}
\usage{
wgtMatrix(g, transpose=TRUE)
adj.check(gm,x)
amat2dag(amat)
find.sink(gm)
}
\details{
  These are not to be called by the user.
  %% or we have been too lazy to document
}
\keyword{ internal }
