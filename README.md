# sCRAP_CommandLine
<b>selective Cross Reactive Antigen Presentation (command line version)</b> 


Commandline version of sCRAP
<br>
https://marisshiny.research.chop.edu/sCRAP/

<b>Please cite this article:</b>
<br>
Yarmarkovich, M., Marshall, Q.F., Warrington, J.M. et al. Cross-HLA targeting of intracellular oncoproteins with peptide-centric CARs. Nature (2021). https://doi.org/10.1038/s41586-021-04061-6

<b>Clone or download the entire repository. Executed the application with R on unix/linux.</b>

```
Usage: Rscript sCRAP.R -p <peptide> -l <HLA-allele> {arguments}

Options:
	-p PEPTIDE, --peptide=PEPTIDE
		(Required) peptide query - {character}

	-l HLA, --hla=HLA
		(Required) hla allele in format: HLA-{allele},  eg.'HLA-A0201' 'HLA-A0301' 

	-s HOTSPOTS, --hotspots=HOTSPOTS
		hotspots locations, numeric value 1-9 seperated by commas

	-o OUTDIR, --outdir=OUTDIR
		Output directory name

	-f FILENAME, --filename=FILENAME
		Output file basename

	-d DATADIR, --datadir=DATADIR
		data directory

	-m HLADIR, --hladir=HLADIR
		hla data directory

	-t FILETYPE, --filetype=FILETYPE
		select output file format (csv/tsv)

	-h, --help
		Show this help message and exit


