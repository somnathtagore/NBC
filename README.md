# NBC Algorithm for Oncogene/Tumor Suppressor Classification

In this study, we use proteins network, mutation and methylation data of fusions and devised a novel network-based parameter called ‘preferential attachment score’ to categorize genes into oncogenes or tumor suppressors. The classification is done using a Naïve Bayes Computation approach. We use an ABC-MCMC method for selecting features for training our classification algorithm. We have performed a survey of tumor suppressors and oncogenes from the perspectives of somatic mutations and network properties. For comparative purposes, we choose some currently well-established methods and found that our algorithm’s efficiency outperforms them with 93.3%. 

There are two scripts, namely, "feature-filter.R" and "NBC.R". 

1. The first script, "feature-filter.R" can be used to test whether the features that need to be used for training are well-defined, feasible and quality-wise good or not. For this purpose, a test file, "feature_filter_data.csv" is provided. 

2. Similarly, the second script "NBC.R" can be used to cluster the genes as OG, TS or DR/PA. For this a test file, "NBC_test_data.csv" is provided.
