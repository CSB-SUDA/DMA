# DMA
Identifying Lymph Node Metastasis-related Factors in Breast Cancer using Differential Modular and Mutational Structural Analysis.
## 
This code repository is all the R language scripts used for analysis.

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
##
##
The video at this address for detailed code execution:http://121.5.52.202/MDAvideo/.

**The R scripts is developed by Xingyi Liu, for questions and comments please contact email--liuxingyisw@outlook.com**
