/*
 * Simple header file including armadillo headers AND causing armadillo to
 * use BLAS and LAPACK by setting the appropriate compiler macros
 *
 * @author Alain Hauser
 * $Id: armaLapack.hpp 13 2011-10-31 16:49:25Z alhauser $
 */

#ifndef ARMALAPACK_HPP_
#define ARMALAPACK_HPP_

#define ARMA_USE_BLAS
#define ARMA_USE_LAPACK
//#define ARMA_DONT_PRINT_RUNTIME_ERRORS
#include <armadillo>

#endif /* ARMALAPACK_HPP_ */
