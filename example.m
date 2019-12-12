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


folder = './data/synthetic/';

addpath(genpath('./accv2016_ba_code'));
addpath(genpath('./iccv2015_sfm_code_mwe'));
addpath(genpath('./sfm'));

imageFiles =   {'lf_capture_17_17_1.mat',...
                'lf_capture_17_17_2.mat',...
                'lf_capture_17_17_3.mat',...
                'lf_capture_17_17_4.mat',...
                'lf_capture_17_17_5.mat',...
                'lf_capture_17_17_6.mat',...
                'lf_capture_17_17_7.mat',...
                'lf_capture_17_17_8.mat',...
                'lf_capture_17_17_9.mat',...
                'lf_capture_17_17_10.mat',...
                'lf_capture_17_17_11.mat',...
                'lf_capture_17_17_12.mat',...
                'lf_capture_17_17_13.mat',...
                'lf_capture_17_17_14.mat',...
                'lf_capture_17_17_15.mat',...
                'lf_capture_17_17_16.mat',...
                'lf_capture_17_17_17.mat',...
                'lf_capture_17_17_18.mat',...
                'lf_capture_17_17_19.mat',...
                'lf_capture_17_17_20.mat',...
                'lf_capture_17_17_21.mat',...
                'lf_capture_17_17_22.mat',...
                'lf_capture_17_17_23.mat',...
                'lf_capture_17_17_24.mat'}

LF = cell(length(imageFiles),1);
H = nan(5,5,length(imageFiles));
D = cell(length(imageFiles),1);

for j=1:length(imageFiles)
    LF{j} = load([folder imageFiles{j}]);
    H(:,:,j) = LF{j}.LF.H;
    D{j} = LF{j}.LF.disparity;
    LF{j} = uint8(255 * LF{j}.LF.LF);
    
    fprintf('%s\n', imageFiles{j})
end
%% Parameters 
inlier_threshold_M = 10e-6;
inlier_threshold_SFM = (10e-6);
p = 0.9999;
ba_iter = 10;
f = 1; 
%% Structure from Motion
% for more details about this part check out our work
% O. Johannsen, A. Sulc, B. Goldluecke : On Linear Structure from Motion for Light Field Cameras
[matches] = lf_extract_matches(LF,H,...
    'inlier_threshold', inlier_threshold_M,...
    'p',p,...
    'support',9,...
    'cross',true);

[pose_A0,matches_inlier] = lf_sfm_sequence_ransac(matches,f,H,...
    'inlier_threshold',inlier_threshold_SFM,...
    'p',p);

%% Bundle Adjustment 
% for more details about this part check out our work
% O. Johannsen, A. Sulc, N. Marniok, B. Goldluecke : Layered scene reconstruction from multiple light field camera views
%this make take a while... around 1.5h
[X_0,rays] = lf_prepare_rays_for_ba(matches_inlier,pose_A0,f);
[X_n, pose_An,f_n] = lf_ba(X_0,pose_A0,rays,f,'niter',25);

%% Unprojecting depth map with 
figure;
axis equal
axis vis3d
hold on;
plot3(X_0(1,:),X_0(2,:),X_0(3,:),'bo')
plot3(X_n(1,:),X_n(2,:),X_n(3,:),'ro')

colors='rgbcmyrgbcmyrgbcmyrgbcmy'; %so beautiful ;)
for i=1:2%length(imageFiles)
    [X,Y,Z] = lf_unproject_dmap(D{i},H(:,:,i),f,[8,8],'step',10);
    [Ri,ti] = invert_Rt(pose_An(1:3,1:3,i),pose_An(1:3,4,i));;
    P = [Ri,ti] * hom([X(:)'; Y(:)';Z(:)']);
    plot3(P(1,:),P(2,:),P(3,:),[colors(i) '.'])
end
