\name{pcalg-internal}
%
%\alias{checkEdges} not exported
\alias{adj.check}
\alias{amat2dag}
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
\alias{lm.cov}
\alias{minUncovCircPath}
\alias{minUncovPdPath}
\alias{my.SpecialDag}
\alias{new.ord}
\alias{pc.cons.intern}%-. ./pc.cons.intern.old
\alias{pcAlgo.Perfect}
\alias{pdsep}%-> ./pdsep.old
\alias{possibleDe}
\alias{qreach}%-> ./qreach.old
\alias{reach}
\alias{rfci.vstructures}
\alias{skeleton.dag2pag}
\alias{triple2numb}
\alias{udag2apag}%-> ./udag2apag.old
\alias{udag2pag} %-> ./udag2pag.old
\alias{udag2pdagRelaxed}%-> ./udag2pdagRelaxed.old
\alias{udag2pdagSpecial}%-> ./udag2pdagSpecial.old
\alias{updateList}
\alias{visibleEdge}
\alias{wgtMatrix}
% FIXME: *do* document [says Martin]  (((or confirm current status))):
\alias{gAlgo-class}
\alias{trueCov}% used in examples / tests!
\alias{allDags}%<- mentioned in other help pages !!
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
