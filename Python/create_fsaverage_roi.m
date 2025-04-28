% Builds a freesurfer overlay .mgh file to define a
% cortical region of interest (ROI).
%
% This follows a two-step procedure:
% 1. Define a set of BA labs
% 2. Intesect with a set of manually defined vertices
% 3. Write .mgh file
%
% This script is part of the repository for our protocol 'Precise
% motor-mapping with transcranial magnetic stimulation'.
% Feel free to adjust to your needs.
%
% Numssen, O., Weise, K., Kalloch, B., Zier, A. L., Thielscher, A.,
% Hartwigsen, G., & Knösche, T. R. (2022, June 9).
% Precise motor-mapping with transcrianial magnetic stimulation - data and code.
% https://doi.org/10.17605/OSF.IO/MYRQN

addpath(genpath('//wsl.localhost/Ubuntu-22.04/usr/local/freesurfer/7.4.1/matlab'))

fs_subject_dir = '//wsl.localhost/Ubuntu-22.04/usr/local/freesurfer/7.4.1/subjects'; %getenv('SUBJECTS_DIR');
side = 'left';

roi_fn = sprintf('%shandknob_.overlay', side);
%% Step 1 -- Brodmann Areals
% Which BA areas to use
%areas = {'Brodmann.6', 'Brodmann.4', 'Brodmann.3', 'Brodmann.1', 'Brodmann.2', 'Brodmann.8', 'Brodmann.44' , 'Brodmann.9'};
areas = {'Brodmann.6', 'Brodmann.4', 'Brodmann.3', 'Brodmann.1', 'Brodmann.2', 'Brodmann.5'};
[vertices, label, colortable] = read_annotation([fs_subject_dir '/fsaverage/label/lh.PALS_B12_Brodmann.annot']);

idx_areas=zeros(size(areas));
for i=1:length(areas)
    for j=1:length(colortable.struct_names)
        if strcmpi(areas{i},colortable.struct_names{j})
            idx_areas(i)=j;
            break
        end
    end
end

% Get indices of BAs
overlay = zeros(size(vertices));

ourvertices = false(size(vertices));
idx_label = colortable.table(idx_areas,5);
for i=1:length(idx_label)
    ourvertices(label==idx_label(i))=true;
end
ourvertices=vertices(ourvertices)+1; % matlab indexing

overlay(ourvertices)=1;

% Step 2 -- Manual reference points
% cut out area around handknob
max_dist = 60; % maximal distance in [mm] to one of the reference points
anti_dist = 130;

% Provide vertex coordinates in coordinate space of the correct file (pial,
% inf, etc).
refpts = [ 12 -4 69; -8 -22 56; -25 -4 39; -11 19 54]; % For left
antirefpts = [100 -10 -40];
% For right:
if strcmpi(side, 'right')
    refpts = [-1 1 1] .* refpts;
    antirefpts = [-1 1 1] .* antirefpts;
end

% Alternatively, you can provide the vertex indices, for example from
% manual selection from FreeView, and grab the coordinates like this:
% vertex_coords(idx,:)

if strcmpi(side, 'left')
    [vertex_coords, faces] = read_surf([fs_subject_dir '/fsaverage/surf/lh.inflated']);
else
    [vertex_coords, faces] = read_surf([fs_subject_dir '/fsaverage/surf/rh.inflated']);
end

allowedVertices=false(size(vertices));
for i=1:size(refpts,1)
    allowedVertices = allowedVertices | sum(bsxfun(@minus, vertex_coords, refpts(i,:)).^2, 2) <= max_dist^2;
end


forbiddenVertices = false(size(vertices));
for i=1:size(antirefpts,1)
    forbiddenVertices = forbiddenVertices | sum(bsxfun(@minus, vertex_coords, antirefpts(i,:)).^2, 2) <= anti_dist^2;
end

% Intersection of BA and reference points
overlay=overlay & allowedVertices & ~forbiddenVertices;
disp([num2str(sum(overlay)) ' vertices found in ROI.' ]);

%%
figure;
ts = trisurf(faces+1, vertex_coords(:,1), vertex_coords(:,2), vertex_coords(:,3)); 
ts.EdgeColor='none'; 
ts.FaceColor = [0.95 0.95 0.95];
hold on

face_overlay = all(ismember(faces+1, find(overlay)), 2);

ts2 = trisurf(faces(face_overlay,:)+1, vertex_coords(:,1), vertex_coords(:,2), vertex_coords(:,3)); 
ts2.EdgeColor='none'; 
ts2.FaceColor = [0.95 0.85 0.5];

xlabel('x'); ylabel('y'); zlabel('z')
axis equal
axis vis3d

light(Style="infinite",Position=[-10 0 10],Color="white");
light(Style="infinite",Position=[0 0 -10],Color=[0.1 0.4 0.8]);


%% Write file
disp(['Writing to ' pwd '/' roi_fn]);
write_curv(roi_fn, overlay, length(vertices));



