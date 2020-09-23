#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::NumericMatrix specDistMatrix(Rcpp::List specList, Rcpp::Function cb)
{
    Rcpp::NumericMatrix ret(specList.length(), specList.length());
    const size_t len = specList.length();
    
    for (size_t i=0; i<len; ++i)
    {
        for (size_t j=i+1; j<=len; ++j)
        {
            // // rows we will operate on
            // Rcpp::NumericMatrix::Row row1 = mat.row(i);
            // Rcpp::NumericMatrix::Row row2 = mat.row(j);
            // 
            // // compute the average using std::tranform from the STL
            // std::vector<double> avg(row1.size());
            // std::transform(row1.begin(), row1.end(), // input range 1
            //                row2.begin(),             // input range 2
            //                avg.begin(),              // output range 
            //                average);                 // function to apply
            // 
            // // calculate divergences
            // double d1 = kl_divergence(row1.begin(), row1.end(), avg.begin());
            // double d2 = kl_divergence(row2.begin(), row2.end(), avg.begin());
            
            // write to output matrix
            double d = Rcpp::as<double>(cb(specList[i], specList[j-1]));
            ret(i, j-1) = d;
            ret(j-1, i) = d;
        }
    }
    
    return ret;
}
