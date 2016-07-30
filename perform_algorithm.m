close all;
clear all;
clc;

%select the correct path: it will open the images' folder
path = '~/Desktop/';

%allows the user to select the first image to analyze
[FileName,PathName] = uigetfile('*.TIF','Select the first image',path);
img1 = sprintf('%s%s',path,FileName);

%allows the user to select the second image to analyze
[FileName,PathName] = uigetfile('*.TIF','Select the second image',path);
img2 = sprintf('%s%s',path,FileName);

%change the value of the following parameters if the user wants different compression
%matrices or metric:
% possible values: tc -> 'Gauss'/'Bernoulli'
%                : f_metric -> 'SIFT'/'SURF'
tc = 'Gauss';
f_metric = 'SIFT';

%run the main function.
main(img1,img2,tc,f_metric);