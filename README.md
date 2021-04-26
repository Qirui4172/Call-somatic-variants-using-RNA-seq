# Call-somatic-variants-using-RNA-seq

This is a pipeline for calling somatic variants using RNA-seq.

The input samples include tumor samples 1, 2, 3 and paired normal samples from mouse. Reference genome mm10.fa and annotation file mm10.gtf were used.

The variants called from RNA-seq contain many false positives, so a further filtration is needed. The suggested filtrations include: read depth in both tumor and normal samples >= 10; MAF > 0.1 in tumor, MAF < 0.01 or 0.02 in normal tissue; strand filtration; base quality & mapping quality > Q20.
