# Molecular Clock Phylodynamics Project
This is the documentation of phylodynamic and transmission dynamics of COVID-19 XBB 1.5 in Indonesia. 

## Data
The data were gathered from GISAID. where we filtered variants 'XBB 1.5' and between Oct 2022-Dec 2023. For outside Indonesia, we take the first cases of XBB 1.5. The reference was taken from GenBank, Wuhan Sequence.

## Method

### Alignment
We using MAFFT in Ubuntu with Linux

```{bash}
mafft -s sample.fasta > sample_alignment.fasta 
```

### Bootstrapping Tree
Ultrafast bootstrapping in IQ-TREE is a technique used to assess the reliability of individual branches in a phylogenetic tree, offering faster computations than traditional bootstrapping methods. By using the -bb 1000 option, it performs 1000 bootstrap replicates to provide statistical support for the tree's branches, enhancing both efficiency and accuracy.

```{bash}
iqtree -s sample_alignment.fasta -m MFP -bb 1000
```
