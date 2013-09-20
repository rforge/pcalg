\name{pcalg-internal}
%
\alias{adj.check}
\alias{allDags}
\alias{amat2dag}
\alias{ci.test}
\alias{find.sink}
\alias{flipEdges}
\alias{gAlgo-class}%-> *do* document! (FIXME)
\alias{labelEdges}
\alias{legal.psep}
\alias{lm.cov}
\alias{make.edge.df}
\alias{orderEdges}
\alias{pcAlgo.Perfect}
\alias{reach}
\alias{trueCov}
\alias{wgtMatrix}
\alias{check.new.coll}
\alias{checkTriple}
\alias{faith.check}
\alias{triple2numb}
\alias{updateList}
\alias{min.uncov.pd.path}
\alias{find.min.discr.path}
\alias{min.uncov.circ.path}
\alias{find.unsh.triple}
\alias{rfci.vstructures}
\alias{dep.triple}
\alias{CheckEdges}
\alias{ancTS}
\alias{skeleton.dag2pag}
\alias{new.ord}
\alias{udag2pdagRelaxed}
\alias{udag2pdagSpecial}
\alias{pc.cons.intern}
\alias{pdsep}
\alias{dreach}
\alias{qreach}
\alias{udag2pag}
\alias{udag2apag}
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
labelEdges(amat)
make.edge.df(amat)
orderEdges(amat)
}
\details{
  These are not to be called by the user.
}
\keyword{ internal }
