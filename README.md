# folks
The folded k-spectrum kernel: a machine learning approach to detecting transcription factor binding sites with gapped nucleotide dependencies


MAIN PROGRAM

The main program is implemented in Matlab, pipelineSVMFE.m. For a given set of input sequences it executes the proposed SVM method in the accompanying paper. The program mainly depends on four Matlab functions which are self-contained and documented, mappingFunction.m, mappingFoldedKspec.m, mappingKspec.m, and binarySpace.m. It also depends on the Tomtom program where the instructions can be found at http://meme-suite.org/tools/tomtom. The executable program tomtom should be placed within the same directory.

Inputs:

• positive sequences, a set of DNA sequences given in fasta format

• r cutoff, a threshold value to filter out weak feature enrichments, 0.005 (default)

Outputs:

• features final, a set of enriched gapped k-mer features 

• ri final, the corresponding enrichment scores, ri


-- For more detail, please see Readme.pdf --
