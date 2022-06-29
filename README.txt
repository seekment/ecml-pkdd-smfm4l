Package: M4Lsmf
Type: Package
Title: Do Small Loss must  Give Large Weight? Avoid Sparse Data Interference with Selective Matrix Factorization
Version: 1.0
Description: This package implements the M4Lsmf algorithm with a matrix tri-factorization framework.
        
Depends:
    MATLAB (>= 2012a)
License: All source code is copyright, under the Artistic-2.0 License.
		For more information on Artistic-2.0 License see [http://opensource.org/licenses/Artistic-2.0](http://opensource.org/licenses/Artistic-2.0)

Files:

M4Lsmf_demo.m: The main function.

M4Lsmf.m: The algorithm for M4Lsmf.

NNDSVD.m: Init the basis matrix factor G with SVD

getOptimalWrWeightsCalculate the weights for inter-associations according to the reconstruction loss and alpha.

getOptimalWhWeightsCalculate the weights for intra-associations according to the trace operation and beta.

