clear raw img img2 img3 img4 block M H H_avg tmp4 xx yy zz img_tmp img_threshold_smooth
save([char(datetime('now','TimeZone','local','Format','yy.MM.dd-HH.mm')) '-Workspace.mat']);  % Dataset images are too heavy to be saved, so reload the raw data and compute them again when opening a workspace.
% Only needed if wanna keep on working, to reopen the raw image that was deleted:
TFM_1_open_stack
