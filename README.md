# LIMITS
# Language: R
# Input: CSV (normalized abundances)
# Output: CSV (correlations)
# Tested with: PluMA 1.0, R 3.2.5

PluMA plugin that identifies important taxa using the LIMITS algorithm
(Fisher and Mehta, 2014) which can infer correlations and magnitudes
using time-series data.

Input CSV file contains a matrix of samples and taxa abundances,
where entry (i, j) is the abundance of taxon j at timepoint i.

The output correlations and values is sent to a CSV file also
as a matrix, where entry (i, j) is the inferred correlation between
taxon i and taxon j.

The plugin source code includes open-source code from the SeqTime package,
which is available under BSD open source license here:
http://hallucigenia-sparsa.github.io/seqtime/

A copy of this BSD license has been included, which covers everything outside of the input(), run() and output()
standard functions for PluMA plugins.  

The professional citation
for SeqTime is here:
Karoline Faust, Franziska Bauchinger, Béatrice Laroche, Sophie de Buyl, Leo Lahti, Alex D Washburne, Didier Gonze, Stefanie Widder, “Signatures of ecological processes in microbial community time series”, Microbiome (2018). 
