clear
clc

% Load Simulated Dataset
save_filename = 'RightBreastMRI.mat';
siminfo_filename = ['sim_info/', save_filename];
load(siminfo_filename,'numRows');
load(siminfo_filename,'numElemPerRow');
numElements = numRows*numElemPerRow;

% Launch Simulations for Cylindrical Wave Transmits
startjob = 1;
for tx_elmt_idx = startjob:numElemPerRow
    GenRFDataSingleCylindricalTxKWave(siminfo_filename, tx_elmt_idx);
end

% Cannot actually assemble all data -- too large
% Must Carefully Choose Desired Frequencies
df = 0.05e6; fmin = 0.2e6; fmax = 0.8e6;
fDATA = fmin:df:fmax;

% Assemble and Save Data from Simulations at Desired Frequencies
AssembleSimData(save_filename, fDATA)