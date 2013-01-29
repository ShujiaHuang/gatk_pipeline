This app applies a number of recalculations to a
[Mappings](http://wiki.dnanexus.com/Types/Mappings) object in order to improve
the quality of variant calls. The pipeline follows the Broad Institute's
recommendations for best practices in variant calling.

Mappings are deduplicated using Picard, realigned around sites of known indels,
and their quality is recalibrated by looking at covariance in quality metrics
with frequently observed variation in the genome. After recalibration, variants
are called with the GATK UnifiedGenotyper module.
