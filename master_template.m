% Location of the scripts:
addpath(genpath('copy_TFM_release1_location_here'));
addpath(genpath('copy_BioFormats_location_here'));

% Information about the video to be analyzed:
% (Prepare the video including the reference frame in fiji, 
% exported as a .tif 8 bit composite hyperstack) 
input_stack='my_raw_data.tif';
n_colors=2;         % Total number of color channels
n_timepoints=2;     % Total number of timepoints including the reference frame
tracking_channel=1; % The channel which contains the matrix tag (e.g. fluorescent nanoparticles)
channel_cell=2;     % The channel which contains cell data (e.g. membrane targeted GFP or cell tracker)
timestep=2;         % Only relevant with time-lapse data
timeunit='min';     % Only relevant with time-lapse data
pxsize_xy=0.2;      % in um (micrometers)
pxsize_z=0.4;       % in um

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell segmentation parameters:
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
% Parameters for the matrix displacement tracking by cross-correlation and for the filtering: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
box_size_xy=41;     %px, odd number.
box_size_z=41;       %px, odd number.
step_size_xy=21;
step_size_z=21;
max_distance_correlation_xy=20;  % px
max_distance_correlation_z=20;   % px
reference_frame=2;%floor(n_timepoints/2);  % Number of the reference frame, or -1 for relative to previous timepoint. 
% Cross-correlation is between -1 (opposite) through 0 (random, not
% matching) to 1 (perfect matching). This parameter sets the cross-correlation
% under which results are rejected for lack of confidence in the matching:
autocorrelation_threshold=0.7;
% Subtract median of the layer to each layer (to compensate for drif):
subtract_median_layer=1;
% In the time direction. Should compute integrals typically if the reference frame is the previous frame.
compute_integrals=0;    
% If computing integrals or if the reference frame is not a relaxed frame,
% should subtract the time-median position considered as the relaxed one:
subtract_median_position=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for the display of the experimental displacement: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% scale factor for displacement vectors display:
f_U=2;      % Same option used for the fitted field. 
% Scale factor picture (if it becomes too big to fit on the screen, creates bugs):
f_pic=3/4;    % Same option used for the fitted field. 
% Should the overlay file be saved on harddrive:
save_overlay_experimental=1;
show_movie_experimental=0;
single_layer=0;  % Same option used for the fitted field. 
layer_saved=147; % Same option used for the fitted field. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for the displacement fitting by force-reconstruction: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Timepoints to process (default: 1:n_timepoints but can be reduced to test
% parameters fast on just one or a few timepoints):
timepoints_to_process=1;%:n_timepoints;%[1 15 30];
% Weight parameters for the minimization problem:
alpha_z=1;                  % Relative importance of good fitting on z values (value between 0->1 if z values should be ignored -> considered as reliable as xy values)
alpha_F_cell_surface=1e-9;  % Relative importance of minimizing the forces on the cell surface 
alpha_F_matrix=1e9;         % Relative importance of minimizing the forces compared to fitting the data in the matrix
alpha_smoothing=5;
% Mechanical parameters:
G=70;      % Shear modulus, Pa.
nu=0;       % Poisson ratio

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Parameters for the display of the reconstructed fields:           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% scale factor for force vectors
f_F=10;
% Should the movie be shown in matlab
show_movie_fitted=0;
% Should the file be saved on harddrive (will be same layer as previously defined for the experimental field):
save_overlay_fitted=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Parameters for the 4D display of the reconstructed fields          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_fraction=0.05;    % fractions of the vectors to be plotted. Start from low values or risk of crash. 
patch_reduction=0.01; % fraction of the facets to be kept in the meshing of the cell surface (initial density is close to the pixel density). 
f_U_3D=3;            % scale factor on the vector display. 
gaussian_smoothing_stress_display=26;    % in um. Gaussian smoothing is applied to the field of tensors on the surface of the cell before using them to find the normal component. 
type_of_color_display=1;    % 0: Norm of the pressure applied by the cell |stress*normal|. 1: Normal pressure normal*stress*normal 2: Principal stress (max |eigenvalue| of stress tensor). 
normalize_color=1;    % If set to 1: normalize to the value below. If set to 0: normalize to max. 
inverted_color_display=0;   % Default (0) is blue for low values, yellow for high values. Set to 1 to invert this. 
minimum_color=-40;
maximum_color=40;          % Value of the contraction used for the max color. % A good guess for this max color is mean(principal_stress(cell_surface),'omitnan')+std(")
fr=4;   % frame rate (not an input, just used below)
tp=180;  % frames per rotation (not an input, just used below)
nt=n_timepoints; % (not an input, just used below)
% The time and rotation is defined for any number of keyframes. The frames in between take a linearly interpolated time and view.
% % % Video from front view and no rotation nor interpolation:
% key_frames=    [1  nt]; 
% key_times=     [1  nt];
% key_azimuths=  [40 54];
% key_elevations=[70 50];
% % Three repeats during full rotation:
% key_frames=    [1  nt*fr  nt*fr+1  2*nt*fr 2*nt*fr+1 3*nt*fr]; 
% key_times=     [1  nt     1       nt      1         nt     ];
% key_azimuths=  [46 46+120 46+120+1 46+240  46+240+1  46+360 ];
% key_elevations=[52 65     65       42      42        52     ];
% % Just one timepoint:
key_frames=1;
key_times=1;
key_azimuths=46;
key_elevations=-52;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for the layer by layer display of the mesh dependent fields  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display of the local pressure and normals are examples of such mesh
% dependent fields
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
pressure_definition=type_of_color_display;    % 0: Norm of the pressure applied by the cell |stress*normal|. 1: Normal pressure normal*stress*normal 2: Principal stress (max |eigenvalue| of stress tensor). 
user_box=[30 40 10 25 15 30];   %xmin xmax ymin ymax zmin zmax
TFM_6b_extract_pressure_data_in_box;
P_box
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             To save the current figure at high resolution:              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%print(['output3D_tif\ManualExport.tif'],'-dtiff','-r600');
