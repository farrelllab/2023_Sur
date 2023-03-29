# Persistent transcriptional states

The scripts enlisted in this folder is to calculate persistence of transcriptional states during zebrafish development, i.e., how long transcriptionally similar cells can be found during development. We used a constant distance ('epsilon') in gene expression space to represent transcriptional similarity. For each cell, cells within that distance was identified and the mean of absolute stage difference between the cells were computed. This code is to optimize the constant distance that was used for the analysis. The underlying principle behind this code is represented as a diagram in Figure 1E.

![Schematic showing our approach for identifying transcriptionally similar cells using an epsilon neighborhood approach and determining whether each cell was in a ‘short-term’ or ‘long-term’ state (based on the mean of absolute stage difference between the analyzed cell and its epsilon-neighbors)](./03-Persistent/persistent_states.png)
