# EvoluSec: Evolutionary Analysis of Protein Secondary Structures

This repository contains the software associated with the publication:  
**Evolutionary model of protein secondary structure capable of revealing new biological relationships (publicly introduced in ISMB 2017, publicly released in February 2019 and formally published in May 2020)**

## Overview

This study pioneered the use of 3D protein structural (variations) evolutinoary model to uncover evolutionary relationships, establishing—for the first time—a functional link between Toll/interleukin-1 receptor (TIR) domains and enzymatic activity

![EvoluSec Diagram](https://github.com/tawssie/EvoluSec/blob/main/image/evolusec.png?raw=true)


The **secondary structure state (DSSP)** evolutionary model was first introduced at the 2017 ISMB conference in the poster:  
**“Inferring protein phylogeny by modelling the evolution of secondary structure”**
The enzymatic activity function was subsequently validated and published in *"NAD+ cleavage activity by animal and plant TIR domains in cell death pathways, Science"*.




## Features

This MATLAB-based software package enables:

- Construction of **phylogenetic trees** based on secondary structure (DSSP)
- Estimation of **evolutionary distances** between proteins based on secondary structure (DSSP)
- Prediction of **ancestral secondary structure states** based on secondary structure (DSSP)

## Purpose

The toolkit is designed to extract evolutionary signals embedded in protein secondary structures—signals that are often overlooked by traditional sequence-based methods. It provides an alternative framework for protein phylogenetic analysis by modelling the evolution of structural features.

## Citation

If you use this software, please cite the following publication:

> Lai, J. S., Rost, B., Kobe, B., & Bodén, M. (2020).  
> *Evolutionary model of protein secondary structure capable of revealing new biological relationships*.  
> *Proteins: Structure, Function, and Bioinformatics*, 88(9), 1251–1259.  
> https://doi.org/10.1002/prot.25898

A publicly available preprint of the Lai et al. manuscript is accessible via bioRxiv:  
> https://www.biorxiv.org/content/10.1101/563452v1

Further enzymatic functions of TIR domains have since been experimentally validated.

> Horsefield, S., Burdett, H., Zhang, X., et al. (2019).  
> *NAD+ cleavage activity by animal and plant TIR domains in cell death pathways*.  
> *Science*, 365(6455), 793–799.  
> https://doi.org/10.1126/science.aax1911

## Dataset

The full dataset associated with the manuscript is publicly available and can be accessed at:  
[Manuscript Dataset](https://drive.google.com/file/d/17f9NTJEef_n6YWleh1jp7Aul5KKja8QP/) (140MB)


## Contributing

Feel free to submit pull requests for improvements.

## History ##
### ISMB 2017 ###
This work established a **protein evolutionary transition probability model** centred on DSSP secondary-structure states, using carefully curated **high-resolution protein crystal structures** for model construction and validation. The methodological framework followed the core logic of classical evolutionary models such as Dayhoff’s PAM and JTT. **Because DSSP states are strongly related to φ/ψ torsion-angle space, the model captures evolutionary transition patterns between secondary-structure states while indirectly reflecting associated changes in local backbone conformation**. Related results were presented as a poster at ISMB 2017 (3DSIG).

<img src="https://github.com/tawssie/EvoluSec/blob/main/image/ISMB_2017_poster.JPG?raw=true" alt="ISMB_poster" width="50%">

### ISMB 2019 ###
Because the number of available high-resolution crystal structures was insufficient to comprehensively explore all joint DSSP–amino acid state transitions, I extended the analysis to a combined structural–sequence representation. At ISMB 2019, I presented a **60-state evolutionary substitution matrix for protein alignment**, defined by three predicted secondary-structure states combined with 20 amino acid states. This work was later incorporated into my PhD thesis, where I further showed that a joint structural-property × amino-acid index can improve protein alignment performance by capturing coupled evolutionary patterns that are not represented by sequence information alone.

<img src="https://github.com/tawssie/EvoluSec/blob/main/image/ISMB_2019_poster.jpg?raw=true" alt="ISMB_poster" width="50%">


**60-state evolutionary substitution matrix alignment and its early NMR-supported SARM1/ARM binding-site context**

My PhD thesis publicly documented an earlier analytical context related to the *Drosophila* SARM1/ARM binding site, integrating the 60-state alignment framework with NMR-related observations.


