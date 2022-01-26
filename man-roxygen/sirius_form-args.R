#' @param relMzDev Maximum relative deviation between the measured and candidate formula \emph{m/z} values (in ppm).
#'   Sets the \option{--ppm-max} command line option.
#' @param elements Elements to be considered for formulae calculation. This will heavily affects the number of
#'   candidates! Always try to work with a minimal set by excluding elements you don't expect. The minimum/maximum
#'   number of elements can also be specified, for example: a value of \code{"C[5]H[10-15]O"} will only consider
#'   formulae with up to five carbon atoms, between ten and fifteen hydrogen atoms and any amount of oxygen atoms. Sets
#'   the \option{--elements} command line  option.
