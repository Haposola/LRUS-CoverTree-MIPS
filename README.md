# LRUS-CoverTree-MIPS
A new tree index structure, named Long Root Unit Sphere Cover Tree (LRUS-CoverTree), and a new k-Maximum Inner-Product Search (k-MIPS) algorithm based on LRUS-CoverTree.

## Dependencies
Depends on MLPACK https://mlpack.org/

## Input data format
any format that compatible with arma::mat.


## command line parameters:
$1 $2 $3 $4 $5 $6 

1: data path.
note: the data folder should have the following three files:
${dataname}_arma.bin, ${dataname}_inside_query.fvecs, ${dataname}_inside_groundtruth.ivecs

2: boolean flag indicating if inside query is used.
inside query: 200 query points are randomly chosen from the original data
random query: 200 query points are randomly generated

3: k, the top-k parameter

4: eps, the approximation quality parameter

5: output filename

6: min_scale of the LRUS-CoverTree

example: E:\\datasets\\audio\\audio 1 50 0.8 E:\\datasets\\audio\\audio_appres.txt -3
in this case, 
the input file: E:\\datasets\\audio\\audio_arma.bin
the query file: E:\\datasets\\audio\\audio_inside_query.fvecs
the ground truth file: E:\\datasets\\audio\\audio_inside_ground_truth.ivecs
