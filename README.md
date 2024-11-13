# Variation in high-amplitude events across the human lifespan

This repository presents a set of code and functions that allows the following by running these codes in the steps:
1. detection of events (get_events.m)
2. randomly sampling matching numbers of events across subjects (with replacement) and detecting clusters (patterns) of events using K-means clustering (randomized_event_sampling_kmeans.m)
3. comparing cluster centroids across iterations to find an optimal cluster centroid that minimizes between-cluster dissimilarity based on the Hungarian algorithm and aligns all partitions to the optimal cluster centroid (get_best_centroid.m).

The test_data.mat file uses publicly available open data from the Midnight Scan Club dataset (Gordon et al., 2017), which is used to identify high-amplitude cofluctations or events in resting state fMRI data. These series of code detects events using edge time series, randomly samples events acriss subjects, creates clusters, aligns clusters to an optimal centroid, and creates a plot of the centroid results.

Clone the repo with this command `git clone https://github.com/joy-neuro/events_lifespan-main.git`

When implementing or referencing this code, please cite the following:
Jo Y., Tanner J., Seguin C., Faskowitz J., Betzel R.F. "Variation in high-amplitude events across the human lifespan." BioRxiv (2024): 2024-05
[https://www.biorxiv.org/content/10.1101/2024.05.15.594378v1.abstract](https://www.biorxiv.org/content/10.1101/2024.05.15.594378v1.abstract)
