# Persistent transcriptional states

The scripts enlisted in this folder is to calculate persistence of transcriptional states during zebrafish development, i.e., how long transcriptionally similar cells can be found during development. We used a constant distance ('epsilon') in gene expression space to represent transcriptional similarity. For each cell, cells within that distance was identified and the mean of absolute stage difference between the cells were computed. This code is to optimize the constant distance that was used for the analysis. The underlying principle behind this code is represented as a diagram in Figure 2A.

![Schematic showing our approach for identifying transcriptionally similar cells using an epsilon neighborhood approach and determining whether each cell was in a ‘short-term’ or ‘long-term’ state (based on the mean of absolute stage difference between the analyzed cell and its epsilon-neighbors)](./persistent_states.png)

To identify the duration of developmental states, we used an epsilon-nearest neighbor approach: we (1) defined a distance in gene expression space to represent transcriptionally similar cells (epsilon), (2) found each cell’s neighbors within that epsilon-neighborhood, and (3) computed the difference in developmental stage between each cell and its neighbors. Cells were categorized based on the mean of absolute stage difference with their neighbors (e.g., <24h, 24–36h, 36–48h, >48h). Cells in transcriptional states that are present “long-term” would be expected to have neighbors that were more different in developmental stage than cells in transcriptional states that occur only during a limited period of development, whose neighbors should all have a similar stage.

In addition, for downstream analyses, cells were classified into “short-term” and “long-term” states based on a threshold of 36 hours. Transcriptional states present ≥36 hours were considered “long-term”. This threshold was chosen to balance focusing on states whose duration was rare while avoiding approaching the upper limit possible for this analysis on a time course of 120 hours. Additionally, each cell was classified as “cycling” or “non-cycling” based on its expression of transcripts associated with different cell cycle phases. 

The code in this folder was used to generate the following figure panels: Figure 2 and Figure S2.

