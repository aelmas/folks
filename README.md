# folks
The folded k-spectrum kernel: a machine learning approach to detecting transcription factor binding sites with gapped nucleotide dependencies

(Please see the file "Readme.pdf" for more detail)

MAIN PROGRAM

The main program is implemented in Matlab, pipelineSVMFE.m. For a given set of input sequences it executes the proposed SVM method in the accompanying paper. The program mainly depends on four Matlab functions which are self-contained and documented, mappingFunction.m, mappingFoldedKspec.m, mappingKspec.m, and binarySpace.m. It also depends on the Tomtom program where the instructions can be found at http://meme-suite.org/tools/tomtom. The executable program tomtom should be placed within the same directory.

Inputs:

• positive sequences, a set of DNA sequences (χi) given in fasta format

• r cutoff, a threshold value to filter out weak feature enrichments, 0.005 (default)

Outputs:

• features final, a set of enriched gapped k-mer features 

• ri final, the corresponding enrichment scores, ri

FEATURE ELIMINATION

For the elimination phase (ri \ rjf ) we considered two options.

Basic: eliminate any feature n with ri(n) > 0, if the corresponding rjf(n) > 0, j = i,...,i+9.

Advanced (default): eliminate any feature n with ri(n) > 0, if it belongs to the "gapped model" of any z with rjf (z) > 0.

In the Advanced option, before the elimination takes place, we run Tomtom for each false enrichment z and retain it if Tomtom finds significant similarity with a particular JASPAR motif (pval < 0.001), otherwise we set rjf(z) = 0, considering them as noise rather than background sequence patterns.
