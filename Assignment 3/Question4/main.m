%% Question 3 implementation
clear;
close all;

%% question a
s50=im2double(imread("slice_50.png"));
s51=im2double(imread("slice_51.png"));
s52=im2double(imread("slice_52.png"));
s50 = padarray(s50, [37, 19], 0);
s51 = padarray(s51, [37, 19], 0);
s52 = padarray(s52, [37, 19], 0);
s50_size = size(s50,1);
s51_size = size(s51,1);
% creating  parallel beam tomographic projections at
% at any 18 randomly angle
rng(4)
r_angles=linspace(0,180,18);
m_slice_50 = radon(s50, r_angles);
m_slice_51 = radon(s51, r_angles);
%%fltered back-projection using the Ram-Lak filter, as implemented in the `iradon' function
reconst_s_50_ramlak= iradon(m_slice_50, r_angles, 'linear', 'Ram-Lak', 1, s50_size);
figure();
imshow(reconst_s_50_ramlak);
colormap('gray')
title("fltered back-projection using the Ram-Lak filter for slice50");


% for slice 51
reconst_s_51_ramlak= iradon(m_slice_51, r_angles, 'linear', 'Ram-Lak', 1, s51_size);
figure();
imshow(reconst_s_51_ramlak);
title("fltered back-projection using the Ram-Lak filter for slice51");
colormap('gray')
%% question b
addpath('./l1_ls_matlab');
questionb(s50,s51);

%% question c
questionc(s50,s51,s52);