%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Location of the scripts                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('copy_TFM_scripts_location_here'));
addpath(genpath('copy_BioFormatsMatlabToolbox_location_here'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Information about the video to be analyzed                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (Video including the reference frame, exported as a .tif 8 bit composite hyperstack from Fiji) 
input_stack='my_raw_data.tif';
n_colors=2;         % Total number of color channels
n_timepoints=2;     % Total number of timepoints including the reference frame
tracking_channel=1; % The channel which contains the matrix tag (e.g. fluorescent nanoparticles)
channel_cell=2;     % The channel which contains cell data (e.g. membrane targeted GFP or cell tracker)
timestep=2;         % Only relevant with time-lapse data
timeunit='min';     % Only relevant with time-lapse data
pxsize_xy=0.2;      % In um (micrometers)
pxsize_z=0.4;       % In um

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Cell segmentation parameters:                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Warning: large values and spherical structuring elements can make the segmentation calculations quickly become untractable. Recommended to start with small values. 
sigma=1;                    % In um, standard deviation of the 3D Gaussian smoothing applied to the image before segmentation.
threshold=0.2;              % Threshold level (ratio of the brightest pixel after smoothing)
fast_segmentation=0;        % If 1, uses a cubic structuring element for the closing and rejection operations (faster, less accurate). If 0, uses a spherical structuring element (slow but accurate). 
closing_radius=2;           % In um. If value above 0, performs a morphological closing operation after thresholding (e.g. to include dark cells in organoids or smooth a cell surface).
fill_holes=1;               % 1 or 0 to activate or not filling of the holes within the thresholded cell(s), e.g. to fill the lumen of organoids or the cytoplasm of a cell imaged with membrane-targeted GFP. 
keep_only_biggest=1;        % Keeps only the largest segmented object, can be useful when there is only one cell/organoid to analyze in the field of view together with small debris to remove.
tolerance_keep_biggest=2;   % In um. Objects in this vicinity from the biggest object are kept together with it (e.g. filopodia or neurites that have small dark areas but still belong to the cell under analysis). 
rejection_radius=2;         % In um. If above 0, experimental displacements closer to the cell(s) than this are automatically rejected, to avoid surface artefacts. If the value of this is the same as the tolerance_keep_biggest, faster computation (done once instead of twice). 
segment_reference_frame=0;  % If the reference frame is a killed cell, no need for segmentation (save time). If the reference frame is just like another frame, and median displacement is then subtracted, then needs to be segmented like any other. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Parameters for the matrix displacement tracking by cross-correlation:  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
box_size_xy=41;                 % in px, odd number. Size of the small 3D boxes tracked in the matrix. 
box_size_z=41;                  % in px, odd number.
step_size_xy=21;                % in px. Spacing of the array on which displacement values are collected. 
step_size_z=21;                 % in px.
max_distance_correlation_xy=20; % in px. Size of the search range (max displacement of a box). 
max_distance_correlation_z=20;  % in px.
reference_frame=2;              % Position of the reference frame, or -1 for tracking relative to the previous timepoint for each frame (not recommended). 
autocorrelation_threshold=0.7;  % Cross-correlation is between -1 (opposite) through 0 (random, not matching) to 1 (perfect matching). This parameter sets the cross-correlation under which results are rejected for lack of confidence in the matching.
subtract_median_layer=1;        % Subtract median displacement of each layer to this layer (to correct for drift).
compute_integrals=0;            % Compute the integral of the displacements in the time direction. Typically used if the reference frame is the previous frame (not recommended).
subtract_median_position=0;     % If computing integrals or if the reference frame is not a relaxed frame, should subtract the time-median position considered as the relaxed one.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Parameters for the display of the overlays:                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_U=2;                          % Scale factor for displacement vectors
f_F=10;                         % Scale factor for force vectors
f_pic=3/4;                      % Scale factor for pictures (if it becomes too big to fit on the screen, creates bugs)
save_overlay_experimental=1;    % Should the experimental displacement overlays be saved on the harddrive
show_movie_experimental=0;      % Should the experimental displacement overlays be shown in Matlab
save_overlay_fitted=1;          % Should the reconstructed fields overlays be saved on the harddrive
show_movie_fitted=0;            % Should the reconstructed fields overlays be shown in Matlab
single_layer=0;                 % Set to 1 to save a single layer of the stack only
layer_saved=0.5;                % Either the number of a layer (integer), or a value between 0 and 1 indicating a position ratio between the bottom and top layer.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Parameters for the displacement fitting by force-reconstruction:     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Timepoints to process (default: 1:n_timepoints-1 to reconstruct forces in
% all frames except the reference one, but can be reduced to test parameters 
% fast on just one or a few timepoints, or to analyze all frames if the median
% displacement is used as a relaxed reference estimate):
timepoints_to_process=1:n_timepoints-1;
G=100;                      % Shear modulus, Pa.
nu=0;                       % Poisson ratio
alpha_F_cell_surface=1e-9;  % Relative importance of minimizing the forces on the cell surface 
alpha_F_matrix=1e9;         % Relative importance of minimizing the forces compared to fitting the data in the matrix
alpha_smoothing=1;          % Relative importance of minimizing the difference between the force in a position of the cell surface and the average force on neighboring cell surface positions
alpha_z=1;                  % Relative importance of good fitting on z values (value between 0->1 if z values should be ignored -> considered as reliable as xy values)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Parameters for the 4D display of the reconstructed fields          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_fraction=0.05;                     % fraction of the vectors to be plotted. Start from low values or risk of crash. 
patch_reduction=0.01;                   % fraction of the facets to be kept in the meshing of the cell surface (initial mesh density is close to the pixel density). 
f_U_3D=3;                               % scale factor on the 3D displacement vector display. 
gaussian_smoothing_stress_display=20;   % in um. Sigma for the Gaussian weighted average applied to the field of tensors on the surface of the cell to interpolate it on every cell surface mesh position. 
type_of_color_display=1;                % 0: Norm of the stress/pressure applied by the cell defined as |stress_tensor*normal|. 1: Normal pressure defined as normal*stress_tensor*normal 2: Principal stress defined as the max eigenvalue of the stress tensor in norm. 
normalize_color=1;                      % If set to 1: normalize to the value below. If set to 0: normalize to max. 
inverted_color_display=0;               % Default (0) is blue for low values, red for high values. Set to 1 to invert this. 
minimum_color=-40;                      % in Pa (same as pN/um^2). Value of the contraction used for the min color. % A good first guess for this min color is min(principal_stress(cell_surface),'omitnan')
maximum_color=40;                       % in Pa (same as pN/um^2). Value of the contraction used for the max color. % A good first guess for this max color is max(principal_stress(cell_surface),'omitnan')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     The time and rotation is defined for any number of keyframes.       %
%          The frames in between are a linear interpolation.              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fr=4;   % frame rate (not an input, just used below)
tp=180;  % frames per rotation (not an input, just used below)
nt=n_timepoints-1; % (not an input, just used below)

% % % Video from front view with no rotation nor interpolation:
% key_frames=    [1  nt]; 
% key_times=     [1  nt];
% key_azimuths=  [50 50];
% key_elevations=[60 60];

% % Three repeats during full rotation:
% key_frames=    [1  nt*fr  nt*fr+1  2*nt*fr 2*nt*fr+1 3*nt*fr]; 
% key_times=     [1  nt     1       nt      1         nt     ];
% key_azimuths=  [46 46+120 46+120+1 46+240  46+240+1  46+360 ];
% key_elevations=[52 65     65       42      42        52     ];

% Just one timepoint:
key_frames=1;
key_times=1;
key_azimuths=46;
key_elevations=-52;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for the layer by layer display of the mesh dependent fields  %
%       (Local pressure and normals are examples of such fields)          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_P=2;
P_fraction=0.2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Running the script execution (comment the parts that you don't want to   %
% execute again when adjusting parameters for the next block to save time) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TFM_1_open_stack;
TFM_2_cell_segmentation;
TFM_3_find_displacement;
TFM_4_displacement_filtering;
%TFM_save_workspace
TFM_4b_display_experimental_displacement_overlay;
TFM_5_force_reconstruction;
plot_force=1;
TFM_5b_display_reconstruction_overlay;
plot_force=0;
TFM_5b_display_reconstruction_overlay;
TFM_6_display3D_and_compute_pressure;
TFM_save_workspace

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      For extraction of the average pressure in a user defined box       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pressure_definition=type_of_color_display;  % Same definitions as above are available.
user_box=[30 40 10 25 15 30];               % xmin xmax ymin ymax zmin zmax. Typically, look at the 3D view, and use the zoom and/or axis() functions to find the coordinates of the box where data should be extracted. 
TFM_6b_extract_pressure_data_in_box;
disp(P_box)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             To save the current figure at high resolution:              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%print(['output3D_tif\ManualExport.tif'],'-dtiff','-r600');
