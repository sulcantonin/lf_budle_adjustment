%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This file is part of the Light Field Structure from Motion package.     %
%                                                                          %
%  This work is licensed under the Creative Commons                        %
%  Attribution-NonCommercial-ShareAlike 4.0 International License.         %
%  To view a copy of this license,                                         %
%  visit http://creativecommons.org/licenses/by-nc-sa/4.0/.                %
%                                                                          %
%  Authors: Ole Johannsen, Antonin Sulc                                    %
%  Contact: first.last@uni-konstanz.de                                     %
%  Website: www.cvia.uni-konstanz.de                                       %
%                                                                          %
%  The package provides code for our ICCV'15 (Structure from Motion) and   %
%  and ACCV'16 paper (adding a bundle adjustment). If you use our work     %
%  we would kindly ask you to cite our below-mentioned papers. Thanks!     %
%                                                                          %
% @string{iccv="IEEE International Conference on Computer Vision (ICCV)"}  %
% @InProceedings{JSG15:iccv,                                               %
%   author = {O. Johannsen and A. Sulc and B. Goldluecke},                 %
%   title = {On Linear Structure from Motion for Light Field Cameras},     %
%   booktitle = iccv, year = {2015}, }                                     %
%                                                                          %
% @string{accv="Asian Conference on Computer Vision (ACCV)"}               %
% @InProceedings{JSMG16:accv,                                              %
%   author = {O. Johannsen and A. Sulc and N. Marniok and B. Goldluecke},  %
%   title = {Layered scene reconstruction from multiple light              %
%       field camera views}, booktitle = accv, year = {2016}, }            %
%                                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

folder = './data/lytro/';

addpath(genpath('./accv2016_ba_code'));
addpath(genpath('./iccv2015_sfm_code_mwe'));
addpath(genpath('./sfm'));

imageFiles =   {'IMG_2292__Decoded.mat',...
                'IMG_2293__Decoded.mat',...
                'IMG_2294__Decoded.mat',...
                'IMG_2298__Decoded.mat'};
       
LF = cell(length(imageFiles),1);
H = nan(5,5,length(imageFiles));
D = cell(length(imageFiles),1);

for j=1:length(imageFiles)
    LF{j} = load([folder imageFiles{j}]);
    H(:,:,j) = LF{j}.RectOptions.RectCamIntrinsicsH;
    LF{j} = uint8(LF{j}.LF(:,:,:,:,1:3) / 255);
    fprintf('%s\n', imageFiles{j})
    
    D{j} = load([folder 'depth/' imageFiles{j}]);
    D{j} = D{j}.disp;
end
%% Parameters 
inlier_threshold_M = 10e-6;
inlier_threshold_SFM = 10e-5;
p = 0.9999;
ba_iter = 10;
%% Structure from Motion
% for more details about this part check out our work
% O. Johannsen, A. Sulc, B. Goldluecke : On Linear Structure from Motion for Light Field Cameras

% find features and their matches
[matches] = lf_extract_matches(LF,H,...
    'inlier_threshold', inlier_threshold_M,...
    'p',p,...
    'support',5,...
    'cross',true);

% estimate f
% this function will plot the residual for different settings of f. Usually a value around 1 works just fine.
lf_estimate_f(matches(1),linspace(0.8,1.6,64));
f=0.94;

[pose_A0,matches_inlier] = lf_sfm_sequence_ransac(matches,f,H,...
    'inlier_threshold',inlier_threshold_SFM,...
    'p',p);

%% Bundle Adjustment 
% for more details about this part check out our work
% O. Johannsen, A. Sulc, N. Marniok, B. Goldluecke : Layered scene reconstruction from multiple light field camera views
[X_0,rays] = lf_prepare_rays_for_ba(matches_inlier,pose_A0,f);
[X_n, pose_An,f_n] = lf_ba(X_0,pose_A0,rays,f,'niter',25);
%% Unprojecting depth map with 
figure;
hold on;
axis equal
axis vis3d
plot3(X_0(1,:),X_0(2,:),X_0(3,:),'bo')
colors='rgbcmyrgbcmyrgbcmyrgbcmy'; %so beautiful ;)
for i=1:length(imageFiles)
    [X,Y,Z] = lf_unproject_dmap(D{i},H(:,:,i),f_n(1),[8,8],'step',10);
    [Ri,ti] = invert_Rt(pose_An(1:3,1:3,i),pose_An(1:3,4,i));;
    P = [Ri,ti] * hom([X(:)'; Y(:)';Z(:)']);
    plot3(P(1,:),P(2,:),P(3,:),[colors(i) '.'])
end
