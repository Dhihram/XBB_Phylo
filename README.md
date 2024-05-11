# Molecular Clock Phylodynamics Project
This is the documentation of phylodynamic and transmission dynamics of COVID-19 XBB 1.5 in Indonesia. The interactive output can be accessed [HERE](https://dhihram.github.io/XBB_Phylo/#Phylogenetic_Tree)

## Data
The data were gathered from GISAID. where we filtered variants 'XBB 1.5' and between Oct 2022-Dec 2023. For outside Indonesia, we take the first cases of XBB 1.5. The reference was taken from GenBank, Wuhan Sequence. In this Github the data are consist:

| Data Name       | Description |
| :--------------| :---------- |
|`metadata_covid.csv`       |metada of samples |
|`sample_alignment.constree`        |the consensus tree|
|`sample_alignment.treefile`           |best-scoring maximum likelihood tree |

## Method

### Alignment
We using MAFFT in Ubuntu with Linux command line

```{bash}
mafft -s sample.fasta > sample_alignment.fasta 
```

### Bootstrapping Tree
Ultrafast bootstrapping in IQ-TREE is a technique used to assess the reliability of individual branches in a phylogenetic tree, offering faster computations than traditional bootstrapping methods. By using the -bb 1000 option, it performs 1000 bootstrap replicates to provide statistical support for the tree's branches, enhancing both efficiency and accuracy.

```{bash}
iqtree -s sample_alignment.fasta -m MFP -bb 1000
```

### Tree dating
Tree dating from tree are using [treedater](https://cran.r-project.org/web/packages/treedater/index.html)

### Transimssion Reconstruction
[epicontact](https://www.repidemicsconsortium.org/epicontacts/) package are using for transmission reconstruction

### Variant Distribution
[outbreakinfo](https://outbreak-info.github.io/R-outbreak-info/) for XBB 1.5 variants distribution
