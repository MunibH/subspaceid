# subspaceid
 
### data
Contains intial data structure. Processed data (`params.doPSTH`) contains `obj.psth`, `obj.condition`, and a separate `meta` struct

### maxdiff_pca 
Subspace ID technique that first identifies 'null' modes via PCA on the preparatory epoch. The projection of the data onto the top null modes are then subtracted from the data. This, ideally, leaves only 'potent' modes within the data. PCA is then used on this data to identify the potent modes

### optimization
Subspace ID technique introduced in [Elsayed 2016](https://www.nature.com/articles/ncomms13239). This method sets up an optimization problem to identify two subspaces (null and potent) such that each subspace maximizes variance within a specified epoch and the two subspaces are orthogonal

### regression
Subspace ID technique introduced in [Kaufman 2016](https://www.nature.com/articles/nn.3643). The method assumed a linear relationship between muscle output and neural activity. The linear operator __W__, once found, contains the null and potent space of neural activity

### utils
Functions used across methods

### getEpochs.m
function to find time indices of specified epochs

### getPsthByCond.m
function to compute psths for each neuron, separated out by specified condition.

### plotLatents.m
function to plot projections of neural data onto a mode

### plotStateSpace.m
function to plot low dimensional trajectories of neural activity

### plotStateSpaceGUI.m 
controls for `plotStateSpace.m`

### subspaceID.m
main function