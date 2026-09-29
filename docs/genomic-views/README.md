# CACTI regional genomic views

Interactive views of the Figure 2B H3K27ac and Figure 2D H3K36me3 examples,
including the focal SNP's +/-500-kb region and three distant control regions.

Coverage is derived from the public datasets analyzed by
[Fair et al. (2024)](https://doi.org/10.1038/s41588-024-01872-x).
See that study's data availability section for the original sequencing and
1000 Genomes genotype sources. Please cite the source study when reusing data.

The exported data include coded sample IDs, focal-SNP genotypes, replicate IDs,
regional coverage, and association annotations. They do not contain full VCFs,
BAM/FASTQ reads, participant names, contact details, or server filesystem paths.
Coded identifiers are linkable to the source datasets, not anonymous labels.

Genotype means weight each plotted profile equally, including technical
replicates. Figure 2B has 78 profiles from 72 donors; Figure 2D has 92 profiles
from 91 donors. Coverage retains the source TMM normalization. BigWigs preserve
base-resolution data; screen rendering can aggregate values when zoomed out.
Regions outside the exports are unavailable, not zero signal.

The viewer and data are static and hosted together. There is no analytics,
form submission, or backend association analysis. Hosting-provider request
logging is separate from the viewer. Reference sequence and external gene
annotations are not loaded.

igv.js 3.0.0 is bundled under its MIT license in `vendor/LICENSE.igv.txt`.
The CACTI code license does not replace terms applicable to source datasets.
