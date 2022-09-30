# Single Cell Analyses of the Mouse Models of Colitis

Analysis of 10x scRNA-seq transcriptomic profile of murine colonic immune cells in control healthy mice, DSS and oxazolone-induced colitis. Our analysis compared the mouse models of colitis to each other, and to a single cell dataset of human colon during health and UC<sup>[1](##References)</sup>, with a focus on myeloid cells.

This repository contains code for:

* Quality control and preprocessing
* MULTI-Seq<sup>[2](##References)</sup> sample demultiplexing using Seurat<sup>[3](##References)</sup>
* Clustering and single cells and marker gene detection
* Ambient RNA contamination removal
* Comparison between control, DSS and oxazolone
* Merging with single cell dataset of human colon during ulcerative colitis
* Comparison between murine and human myeloid cells.

Please contact Ojan Khosravifar if you have any questions (ojankhosrvaifar@gmail.com).

## References
1.   Smillie, C. S. *et al*. Intra- and Inter-cellular Rewiring of the Human Colon during Ulcerative Colitis. *Cell* **178**, 714-730 (2019).
2.   McGinnis, C. S. *et al*. MULTI-seq: sample multiplexing for single-cell RNA sequencing using lipid-tagged indices. *Nature Methods* **16**, 619-626 (2019).
3.   Hao, Y. *et al*. Integrated analysis of multimodal single-cell data. Integrated analysis of multimodal single-cell data. *Cell*, **184**, 3573-3587 (2021).
