%% Get the origin from the NIfTI file
niftiFile = 'C:/Users/bnplab-admin/Documents/MRI/mri_110/NEURO_ZIEMANN_CONNECT2BRAIN_20211103_161821_369000_T1w_MPR_20211103161821_7.nii';

niftiInfo = niftiinfo(niftiFile);
niftiOrigin = niftiInfo.Transform.T(1:3, 4)';

%% Get the origin from the DICOM files
% Specify the path to the folder containing IMA files
dicomFolder = 'B:/MRI Data/DICOM Source Data/TMSEEG_SEG_110_TMSEEG_SEG_110/NEURO_ZIEMANN_CONNECT2BRAIN_20211103_161821_369000/T1W_MPR_0007';

% Get a list of all DICOM files in the folder
dicomFiles = dir(fullfile(dicomFolder, '*.IMA'));

% Preallocate arrays for DICOM origins
dicomOrigins = zeros(length(dicomFiles), 3);

% Loop through each DICOM file to get the origins
for i = 1:length(dicomFiles)
    dicomFilePath = fullfile(dicomFolder, dicomFiles(i).name);
    dicomInfo = dicominfo(dicomFilePath);
    dicomOrigins(i, :) = dicomInfo.ImagePositionPatient';
end

%% Display the origins
disp(['NIfTI Origin: ' num2str(niftiOrigin)]);
% Display the origins of the first DICOM file
disp(['DICOM Origin (First File): ' num2str(dicomOrigins(1, :))]);

% Check if the origins are approximately equal
if all(isequal(round(dicomOrigins, 5), round(dicomOrigins(1, :), 5)))
    disp('Coordinate systems align.');
else
    disp('Coordinate systems do not align.');
end






