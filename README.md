A python script to specify the location of the read from sam file with transcripts information on the mRNA 


mRNA consists of three specific regions 5'UTR CDS 3'UTR


in Mouse genome mm10: 
we found that average length of 5'UTR is 173 nucleotides
we found that average length of CDS is 1128 nucleotides
we found that average length of 3'UTR is 710 nucleotides

so, we created custom number of bins where 5'UTR is represented as 20 bins, CDS is represented as 100 bin, 3'UTR is represented as 70 bins



python align_sam_transcripts_custom_bins.py -SamFile test_genomic.sam -utr3Folder ./utr3/ -utr5Folder ./utr5/ -cdsFolder ./cds/ -output test_custom.tab -sora test_custom.png
