CCC HotSpot Calling
===========



Input Format Requirements
-----------

* VCF
* Sorted
* Indexed


Docker Image Dependencies
-----------

* [Bcftools](https://github.com/ohsu-computational-biology/dockerized-tools/tree/develop/bcftools)
* [VEP](https://github.com/ohsu-computational-biology/dockerized-tools/tree/develop/vep)
* [vcf2maf](https://github.com/ohsu-computational-biology/dockerized-tools/tree/develop/vcf2maf)
* [OncodriveFM](https://github.com/ohsu-computational-biology/dockerized-tools/tree/develop/oncodrivefm)
* [Taylor Lab Algorithm](https://github.com/ohsu-computational-biology/dockerized-tools/tree/develop/taylor-lab-hotspots)

Running Locally
-----------

1. Download [cromwell](https://github.com/broadinstitute/cromwell/releases/download/0.22/cromwell-0.22.jar)
2. Run: `java -jar cromwell-0.22.jar native_cromwell.wdl native_cromwell.json`
