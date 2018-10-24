# folks
The folded k-spectrum kernel: a machine learning approach to detecting transcription factor binding sites with gapped nucleotide dependencies$^1$

MAIN PROGRAM

The main program is implemented in Matlab, pipelineSVMFE.m. For a given set of input sequences it executes the proposed SVM method in the accompanying paper. The program mainly depends on four Matlab functions which are self-contained and documented, mappingFunction.m, mappingFoldedKspec.m, mappingKspec.m, and binarySpace.m. It also depends on the Tomtom program where the instructions can be found at http://meme-suite.org/tools/tomtom. The executable program tomtom should be placed within the same directory.

Inputs:

• positive sequences, a set of DNA sequences given in fasta format

• r cutoff, a threshold value to filter out weak feature enrichments, 0.005 (default)

Outputs:

• features final, a set of enriched gapped k-mer features 

• ri final, the corresponding enrichment scores, ri


-- For more detail, please see Readme.pdf --

$^1$ Elmas A, Wang X, Dresch JM (2017) The folded k-spectrum kernel: A machine learning approach to detecting transcription factor binding sites with gapped nucleotide dependencies. PLoS ONE 12(10): e0185570. https://doi.org/10.1371/journal.pone.0185570
