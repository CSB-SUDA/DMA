# DMA
Identifying Lymph Node Metastasis-related Factors in Breast Cancer using Differential Modular and Mutational Structural Analysis.
## how to use code
This code repository is all the R language scripts used for analysis.
```{r, eval = FALSE}

section-ID. bash fetchType.sh [Required packages] (The input file [from]) (output-file-ID:The output file[annotation])

There requires preparing things, inputfile befor.Then executing "bash fetchType.sh".finally, you will get outputfile.:

Required packages:Executing code requires a pre-prepared environment.
from: NEW/section-ID.output-file-ID NEW:Data that needs to be downloaded, or data that is generated;x.x:The section ID that produces the data.
The input file:Documents that need to be prepared.
The output file:The file that the script outputs.

Files that need to be downloaded:
* gencode.v22.annotation.gtf.gz:https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_22/gencode.v22.chr_patch_hapl_scaff.annotation.gtf.gz
* TCGA-BRCA.survival.tsv:https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-BRCA.survival.tsv
* TCGA-BRCA.htseq_fpkm.tsv.gz:https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-BRCA.htseq_fpkm.tsv.gz
* TCGA-BRCA.GDC_phenotype.tsv.gz:https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-BRCA.GDC_phenotype.tsv.gz
* string_interactions_short_0_4.tsv:https://string-db.org/ <- DEGs from 2.4 combined score > 0.4.
* reactome.homo_sapiens.interactions.psi-mitab.txt:https://reactome.org/download/current/interactors/reactome.homo_sapiens.interactions.psi-mitab.txt.
* ReactomePathways.txt:https://reactome.org/download/current/ReactomePathways.txt
* Cancer_Gene_Census_Hallmarks_Of_Cancer.tsv:https://cancer.sanger.ac.uk/cosmic. login to download.
* TCGA.BRCA.mutect.995c0111-d90b-4140-bee7-3845436c3b42.DR-10.0.somatic.maf.gz:https://portal.gdc.cancer.gov/.

1. bash fetchType.sh 
	@requried:
	@input file:
		gencode.v22.annotation.gtf[NEW]
	@output file:
		1:protein_coding_genecode_v22.csv[]

2. Rscript 1.DifferetialAnalysis.R 
	@requried:
		R:ggplot2 stringr limma ggrepel GeoTcgaData
	@input file:
		./data/TCGA-BRCA.survival.tsv[NEW]
		./data/TCGA-BRCA.htseq_fpkm.tsv.gz[NEW]
		./data/TCGA-BRCA.GDC_phenotype.tsv.gz[NEW]
		./data/protein_coding_genecode_v22.csv[1.1]
	@output file:
		1:BRCA_clinical.csv[]
		2:BRCA_log2TPM_DEG.csv[]
		3:BRCA_log2TPM_DEGmatrix.csv[]
		4:BRCA_log2TPM_matrix.csv[]
		5:volo.pdf[Fig S1a]

3. Rscript 2.survivalAnalysis.R
	@requried:
		R:survival survminer RColorBrewer
	@input file:
		./data/BRCA_clinical.csv[2.1]
	@output file:
		1:./survival.pdf[Table S1b]

4. Rscript 3.NetworkWeighted&modular.R
	@requried:
		R:igraph(1.2.11) reshape2 ggpubr ggplot2 networkD3
		#igraph must be this version.
	@input file:
		./data/BRCA_clinical.csv[2.1]
		./data/BRCA_log2TPM_DEGmatrix.csv[2.3]
		./data/string_interactions_short_0_4.tsv[NEW]
	@output file:
		1:ecN0Clusterlg10.csv[Fig 2a]
		2:ecNplusClusterlg10.csv[Fig S1d]
		3:stringdb.csv[]
		4:sankey.html[Fig 2b]
		5:stringdb.csv[]

5. Rscript 4.Jaccardsimilarity.R
	@requried:
	@input file:
		./data/ecN0Clusterlg10.csv[4.1]
		./data/ecNplusClusterlg10.csv[4.2]
	@output file:
		1:Jaccardsimilarity.pdf[Fig S1c]


6. Rscript 5.edgeInfo.R
	@requried:
		R:survival survminer poolr
	@input file:
		./data/ecN0Clusterlg10.csv[4.1]
		./data/ecNplusClusterlg10.csv[4.2]
		./data/BRCA_log2TPM_DEGmatrix.csv[2.3]
		./data/BRCA_clinical.csv[2.1]
		./stringdb.csv[4.5]
	@output file:
		1:parameters.csv[]

7. Rscript 6.knockDownTest.R
	@requried:
		R:igraph ggplot2
	@input file:
		./data/parameters.csv[6.1]
	@output file:
		1:betweennessN0edge.pdf[Fig 3a]
		2:shortestPathWayN0edge.pdf[Fig 3b]
		3:N0betweennessinterTable.csv[]
		4:N0betweennessintraTable.csv[]
		5:N0closenessinterTable.csv[]
		6:N0closenessintraTable.csv[]


8. Rscript 7.edgeEnrichment.R
	@requried:
		R:biomaRt clusterProfiler
	@input file:
		./data/parameters.csv[6.1]
		./data/reactome.homo_sapiens.interactions.psi-mitab.txt[NEW]
		./data/ReactomePathways.txt[NEW]
	@output file:
		1:interEdge.csv[]
		2:intraEdge.csv[]
		3:interEdge.pdf[Fig 3c]
		4:intraEdge.pdf[Fig 3d]

9. Rscript 8.datehub&partyhub.R
	@requried:
		R:biomaRt clusterProfiler
	@input file:
		./data/parameters.csv[6.1]
		./data/reactome.homo_sapiens.interactions.psi-mitab.txt[NEW]
		./data/ReactomePathways.txt[NEW]
	@output file:
		1:dateHubs.csv[]
		2:partHubs.csv[]
		3:BothChangeEdgDegree.pdf[Fig 4a]

10. Rscript 9.robustnessTest.R
	@requried:
		R:biomaRt clusterProfiler
	@input file:
		./data/parameters.csv[6.1]
		./data/dateHubs.csv[9.2]
		./data/partHubs.csv[9.3]
	@output file:
		1:N0.pdf[Fig 4b]
		2:Nplus.pdf[Fig 4c]

11. Rscript 10.hallmarkEnrichment.R
	@requried:
		R:org.Hs.eg.db clusterProfiler ggplot2
	@input file:
		./data/protein_coding_genecode_v22.csv[1.1]
		./data/dateHubs.csv[9.2]
		./data/partHubs.csv[9.3]
		./data/Cancer_Gene_Census_Hallmarks_Of_Cancer.tsv[NEW]
	@output file:
		1:hallmark.pdf[Fig 4d]
		2:partKEGGData.csv[]
		3:TERM2GENE.csv[]
		4:dateKEGGData.csv[]

12. Rscript 11.dynamicModule.R
	@requried:
		R:igraph GOSemSim tidyverse dplyr EM.R
	@input file:
		./data/protein_coding_genecode_v22.csv[1.1]
		./data/dateHubs.csv[9.2]
		./data/partHubs.csv[9.3]
		./data/Cancer_Gene_Census_Hallmarks_Of_Cancer.tsv[NEW]
	@output file:
		1:TFC_noweight.csv[Fig 5a]

13. Rscript 12.Rosechart.R
	@requried:
		R:ggplot2
	@input file:
		./data/TFC_noweight.csv[12.1]
	@output file:
		1:Rosechart.pdf[Fig 5b]

14. Rscript 13.gsea.R
	@requried:
		R:clusterProfiler enrichplot
	@input file:
		./data/BRCA_log2TPM_DEG.csv[2.2]
		./data/c2.cp.kegg.v7.5.1.symbols.gmt[12.1]
	@output file:
		1:KEGG_CALCIUM_SIGNALING_PATHWAY.pdf[Fig 5c]

15. Rscript 14datahubParthubMutation.R
	 @requried:
		R:maftools ggplot2
	@input file:
		./data/dateHubs.csv[9.2]
		./data/partHubs.csv[9.3]
		./data/BRCA_clinical.csv[2.1]
		./data/TCGA.BRCA.mutect.995c0111-d90b-4140-bee7-3845436c3b42.DR-10.0.somatic.maf.gz[NEW]
	@output file:
		1:datelaml.pdf[Fig 5d]
		2:partlaml.pdf[Fig S3a]

This is the code for the main result of Fig1-4.
Please contact the liuxingyisw@outlook.com if you have any questions.
```
##Brief introduction


In this project, we focused on the rewiring of the co-expression network of differential genes under different conditions (lymph node metastases and groups without lymph node metastases). We compared the network at the level of modules, interactions and nodes, and used date hub and the interactions highlighted by differetial module analysis to identify the core dynamic structure of the network.

### network comparison
First, we constructed two differential gene weighted networks with identical topologies but different weights.
The network weighted process under different conditions is shown in the following figure.

<div align=center>
  <img width="500" src="https://github.com/CSB-SUDA/DMA/blob/main/picture/weightedNetwork.png"/>
</div>
<p align=center>The construction flow of weighted network under different conditions</p>

Due to the strengthening and weakening of the interaction in the network, the module division of the network changes.
<div align=center>
  <img width="800" src="https://github.com/CSB-SUDA/DMA/blob/main/picture/twonetwork.png"/>
</div>
<p align=center>Division of weighted network modules.</p>

What is the biology functional impact of the remodule caused by the rewiring of these modules?We need to annotate these changes at the module level。
<div align=center>
  <img width="400" src="https://github.com/CSB-SUDA/DMA/blob/main/picture/changeDetail.png"/>
</div>
<p align=center>Changes between modules between no-LNM and LNM.</p>

We propose a highlighting method for identifying dynamic core modules using the differential module approach.

<div align=center>
  <img width="1000" src="https://github.com/CSB-SUDA/DMA/blob/main/picture/pipline.png"/>
</div>
<p align=center>Simplified process, detailed process reference article</p>
##


**The R scripts is developed by Xingyi Liu, for questions and comments please contact email--liuxingyisw@outlook.com**
