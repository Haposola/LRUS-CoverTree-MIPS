# LRUS-CoverTree-MIPS
A new tree index structure, named Long Root Unit Sphere Cover Tree (LRUS-CoverTree), and a new k-Maximum Inner-Product Search (k-MIPS) algorithm based on LRUS-CoverTree.

# Dependencies
Depends on MLPACK https://mlpack.org/

# Input data format
any format that compatible with arma::mat.


# command line parameters:
$1 $2 $3 $4 $5 $6
$1:data path. note: the data folder should have the following three files:
${dataname}_arma.bin, ${dataname}_inside_query.fvecs, ${dataname}_inside_groundtruth.ivecs

$2: boolean flag indicating if inside query is used.
inside query: 200 query points are randomly chosen from the original data
random query: 200 query points are randomly generated

$3: k, the top-k parameter

$4: eps, the approximation quality parameter

$5: output filename

%6: min_scale of the LRUS-CoverTree
