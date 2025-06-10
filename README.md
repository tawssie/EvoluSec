# EvoluSec
Evolutionary model of protein secondary structure capable of revealing new biological relationships

├── _README.md          (README for "Evolu-sec-Package")
├── _util/               (Scripts and libraries, categorized by folders)
│   ├── ASSR/            (Ancestral Secondary Structure Reconstruction)
│   │   ├── ASSR.m   (ASSR algorithm)
│   │   ├── CharVector2JavaArray.m
│   │   ├── SetTransitionProb.m
│   │   └── NBbnkit.jar  (Pre-compiled JAVA library to calculate exact maximum likelihood values ASSR)
│   ├── MLDist/          (Evolutionary distance estimation by maximum likelihood)
│   │   ├── AADIST.m
│   │   ├── SSDIST.m
│   │   └── BrentOptAlgo4Dist.p
│   ├── MLTree/          (Tree inference by maximum likelihood)
│   │   ├── CombineOptimization02.m
│   │   ├── CreatePhytreeFromBGraph.m
│   │   ├── GetInitialDistanceByLeastSquare.m
│   │   ├── GetTipsLikelihood4SS.m
│   │   ├── GetTreeLikelihood03.m
│   │   ├── ImproveTree04.m
│   │   ├── NNInterchange03.m
│   │   ├── OptimizeSwap02.m
│   │   ├── SSMLTree.m
│   │   └── BrentOptAlgo4Tree.p
│   └── Model/                   (Amino acid and secondary structure models)
│       ├── DSO.m                (Amino acid model, Dayhoff's)
│       ├── GetExpMByTime.m      (Script to calculate transition probabilities)
│       ├── JTT.m                (Amino acid model, JTT)
│       └── EvoSS.m              (Secondary structure model, DSO-SS, JTT-SS)
├── Evolu-sec-ASSR/              (ASSR scripts with an example)
│   ├── _README.txt
│   ├── DoASSR.m                 (Major script to do ASSR)
│   ├── Main.m
│   ├── example_Wilson_MSSA.fasta      (MSSA from Wilson et al.)
│   └── example_Wilson_PrunedTree.txt  (Phylogenetic from Wilson et al.)
├── Evolu-sec-MLDist/            (Evolutionary distance estimation with examples)
│   ├── _README.txt
│   ├── EstimateAADistance.m     (Script to estimate evolutionary distance based on JTT model)
│   ├── EstimateSSDistance.m     (Script to estimate evolutionary distance based on DSO-SS or JTT-SS)
│   ├── Main.m
│   ├── example_PF00017_AA.fasta (Amino acid example 1)
│   ├── example_PF00017_SS.fasta (Secondary structure example 1)
│   ├── example_PF00031_AA.fasta (Amino acid example 2)
│   └── example_PF00031_SS.fasta (Secondary structure example 2)
├── Evolu-sec-MLTree/            (Evolutionary tree inferences)
│   ├── _README.txt
│   ├── BuildSSTree.m            (Script to construct phylogenetic tree based on Evolu-sec model)
│   ├── Main.m
│   └── example.fasta
├── Evolu-sec-Boundary/          (Variation in boundaries of secondary structure elements)
│   ├── EvoSS_BoundaryTest1.m    (Test boundaries between modern and actual, reconstructed ancestral protein)
│   ├── EvoSS_BoundaryTest2.m    (Test boundaries within phylogenetic trees)
│   └── Main.m
└── _Dataset/                    (All datasets and model parameters)
    ├── ASSR/                    (ASSR-related dataset)
    ├── Evolu-sec_Model          (Evolu-sec model)
    ├── Pfam/                    (MSA, MSSA from Pfam)
    └── TLR1_MLTree/             (TLR1_MLTrees, dataset, GO terms, Dali searching results)

