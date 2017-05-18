# cloningoligos


GOAL: Generate oligos for cloning for many genes in automated fashion using cDNA

INPUT: FASTA file of CCDS sequences for genes of interest, and cDNA FASTA sequences
List of Restriction enzyme sites, in order, in the desired vector of interest.

Align CCDS with associated cDNA. Generate Oligo of at least 20-25 bp at TSS and stop codon. 
Use Primer3 to determine if TM of both oligos is satisfactory: if not, extend the oligo with 
the lowest TM until the TM is within 5 degrees of the second one. 

Scan for possible restriction sites within trasncript. Add compatible sites (not found in sequence) on the oligos.

Add Buffer sequence after RE site of at least 6 bp, to allow direct digestion of PCR product. 


