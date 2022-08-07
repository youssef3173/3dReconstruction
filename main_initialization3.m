clear;
close all;
clc;
addpath(genpath('.'));

%% Load image names 
load('im_names.mat'); %'im_names'

%% Load tracks 
load('tracks.mat');

%% Linear calibration matrix
K = [535 0 320;
    0 539 247;
    0 0 1];

%% Define first pair of images and plot matches 
id_im1 = 1;
id_im2 = 2;
id_im3 = 3;

im1 = double(imread(fullfile('./images/',im_names{id_im1})))/255;
im2 = double(imread(fullfile('./images/',im_names{id_im2})))/255;
im3 = double(imread(fullfile('./images/',im_names{id_im3})))/255;

p1 = p{id_im1};
p2 = p{id_im2};


[Li1, Loc2] = ismember(tracks{id_im1}.point3D_ids, tracks{id_im2}.point3D_ids);
matches = [tracks{id_im1}.point2D_ids(Li1)'; tracks{id_im2}.point2D_ids(Loc2(Li1))'];


%% Initialization
Rw1 = eye(3);
Tw1 = zeros(3,1);

Rw2 = Rw1*(R21');
Tw2 = -Rw2*T21 + Tw1;

[Li1, Loc2] = ismember(tracks{id_im1}.point3D_ids, tracks{id_im2}.point3D_ids);
matches = [tracks{id_im1}.point2D_ids(Li1)'; tracks{id_im2}.point2D_ids(Loc2(Li1))'];

p1 = p{1}(:,matches(1,:));
p1_hom = [p1; ones(1,size(p1,2))];
m1_hom = inv(K)*p1_hom;
p2 = p{2}(:,matches(2,:));
p2_hom = [p2; ones(1,size(p2,2))];
m2_hom = inv(K)*p2_hom;

Uw = triangulate(m1_hom, m2_hom, Rw1, Tw1, Rw2, Tw2);
Rwi_c = {Rw1, Rw2};
Twi_c = {Tw1, Tw2};

im_names_c = {fullfile('./images/',im_names{id_im1}),fullfile('./images/',im_names{id_im2})};

%% Bundle Adjustment

%define tracks currently used and define new point3D_ids
tracks_cur = {tracks{id_im1}, tracks{id_im2}};
nTracks = sum(Li1);
tracks_cur{1}.oldpoint3D_ids = tracks_cur{1}.point3D_ids(Li1);
tracks_cur{1}.point2D_ids = tracks_cur{1}.point2D_ids(Li1);
tracks_cur{1}.point3D_ids = 1:nTracks;


tracks_cur{2}.oldpoint3D_ids = tracks_cur{2}.point3D_ids(Loc2(Li1));
tracks_cur{2}.point2D_ids = tracks_cur{2}.point2D_ids(Loc2(Li1));
tracks_cur{2}.point3D_ids = 1:nTracks;

p1_vis = p{1}(:,tracks_cur{1}.point2D_ids);

maxIt = 25;
[Rwi_BA_c, Twi_BA_c, Uw_BA] = bundleAdjustment(Rwi_c, Twi_c, Uw, tracks_cur, K, p, nTracks, maxIt,2);
%%
% Points déjà construits à partir des images 1 et 2 
trackConstructed = tracks_cur{2}.oldpoint3D_ids'; 
track = tracks_cur{2}.point3D_ids; 
%% Ajout de la caméra 3 
for nCam =3:4
    %% estimation de la position de la caméra 
    [Rwi_BA_c,Twi_BA_c] = getPosiCam(Rwi_BA_c,Twi_BA_c,Uw_BA, tracks_cur, tracks, nCam,K,p); 
    
    % initialisation 
    [Li1, Loc2] = ismember(tracks{nCam}.point3D_ids, trackConstructed);
    tracks_cur{nCam}.oldpoint3D_ids = tracks{nCam}.point3D_ids(Li1);
    tracks_cur{nCam}.point2D_ids = tracks{nCam}.point2D_ids(Li1);
    tracks_cur{nCam}.point3D_ids=track(Loc2(Li1));
    im_names_c = [im_names_c,{fullfile('./images/',im_names{nCam})}];
    % triangulation 
    % on compare à chaque itération im_nCam avec im_i where i in [1:nCam-1] 
    for i=1:nCam-1
        % on ne choisi que les points qui n'existent pas dans 
        % le nuage initial
        [Li_i, ~] = ismember(tracks{i}.point3D_ids, trackConstructed);
        Lia = not(Li_i); 
        % les nouveaux points envisegables 
        new3D = tracks{i}.point3D_ids(Lia); 
        new2D = tracks{i}.point2D_ids(Lia); 
        % L'étape du matching 
        [Li1, Loc2] = ismember(new3D, tracks{nCam}.point3D_ids);
        matches = [new2D(Li1)'; tracks{nCam}.point2D_ids(Loc2(Li1))'] ; 
        
        % On rajoute au tracks_cur les nouveaux points  
        tracks_cur{nCam}.oldpoint3D_ids = [tracks_cur{nCam}.oldpoint3D_ids;tracks{nCam}.point3D_ids(Loc2(Li1))];
        tracks_cur{nCam}.point2D_ids = [tracks_cur{nCam}.point2D_ids;tracks{nCam}.point2D_ids(Loc2(Li1))];
        tracks_cur{nCam}.point3D_ids = [tracks_cur{nCam}.point3D_ids,(length(Uw_BA)+1):length(Uw_BA)+length(new3D(Li1))];
        
        % on réalimente le nuage de points 
        trackConstructed = [trackConstructed,tracks{nCam}.point3D_ids(Loc2(Li1))'];
        track  = [track,(length(Uw_BA)+1):length(Uw_BA)+length(new3D(Li1))];
        
        p1 = p{i}(:,matches(1,:));
        p1_hom = [p1; ones(1,size(p1,2))];
        m1_hom = inv(K)*p1_hom;
        p2 = p{nCam}(:,matches(2,:));
        p2_hom = [p2; ones(1,size(p2,2))];
        m2_hom = inv(K)*p2_hom;

        UwnCam_i = triangulate(m1_hom, m2_hom, Rwi_BA_c{i}, Twi_BA_c{i}, Rwi_BA_c{nCam}, Twi_BA_c{nCam});
        Uw_BA = [Uw_BA, UwnCam_i];
        
    end
end
    
maxIt = 25;
[Rwi_BA_c, Twi_BA_c, Uw_BA] = bundleAdjustment(Rwi_BA_c, Twi_BA_c, Uw_BA, tracks_cur, K, p, size(Uw_BA,2), maxIt,nCam);

%% Visualize refined 3D scene
figure, show3D(Rwi_BA_c, Twi_BA_c, K, Uw_BA, im_names_c); title('3D view after BA');
drawnow;

    
