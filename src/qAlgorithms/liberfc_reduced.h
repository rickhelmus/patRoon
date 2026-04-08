// header to export real-valued functions relevant to qAlgorithms from the liberfc source
// code. This is done to avoid having required dependencies, especially since we currently
// do not operate in the domain of complex numbers

namespace liberfc
{
    double erfi(double x);

    double dawson(double x);
}