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
\alias{faith.check}
\alias{find.min.discr.path}
\alias{find.sink}
\alias{find.unsh.triple}
\alias{flipEdges}
\alias{legal.path}
\alias{minUncovCircPath}
\alias{minUncovPdPath}
\alias{my.SpecialDag}
\alias{pcAlgo.Perfect}
\alias{possibleDe}
\alias{reach}
\alias{rfci.vStruc}% <- renamed; TODO (MM)
\alias{triple2numb}
\alias{updateList}
\alias{visibleEdge}
\alias{wgtMatrix}
%--------------------------------------------------------------------------
% FIXME: *do* document [says Martin]  (((or confirm current status))):
\alias{gAlgo-class}
\alias{trueCov}% used in examples / tests! -- change name !
\alias{allDags}%<- mentioned in other help pages !!
\alias{amat2dag}% in ../tests/test_amat2dag.R
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
