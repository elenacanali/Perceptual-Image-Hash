function [authenticity, blocks_RGB]=main(img1,img2,tc,f_tform)

%% Optional input
% Manages the optional inputs:
% tc : type of compression, default value 'Gauss', alternative 'Bernoulli'
% f_tform : feature extraction, default value 'SIFT', alternative 'SURF'

if (nargin < 3)
   tc = 'Gauss'; 
end
if (nargin < 4)
    f_tform = 'SIFT';
end

%% Load image

image_0=single(im2double(imread(sprintf(img1))));

%% Set thresholds and other parameters

Tr=3.5;     % if max(d)>=Tr then image inauthentic.

Tr_a=0.5;   % if min(d_a)>Tr_a then geometric transformation.
            
Tr_1=4;     % Threshold used in the localization, when no geometric 
            % distortion is found.
Tr_2=4.5;   % Threshold used in the localization, when any geometric 
            % distortion is found.

s1=9.5/100; % projection rate s1
            
s2=9/100;   % projection rate s2

% update thresholds when SURF metric is used
if(strcmp(f_tform,'SURF'))

    Tr_a = 1;
    s1 = 9.5/100;

end

%% Hash generation algorithm
[F_0, Gs1, Gs2]=extract_features(image_0, s1, s2,tc, f_tform);
K=rand(1); % Secret key-seed shared by a sender and a receiver 
F_E=encrypt(F_0, K);

%% Verification algorithm

image_t=single(im2double(imread(sprintf(img2))));
[d, d_a,blocks_RGB]=verifier(image_t, F_E, K, Gs1, Gs2, Tr_a, f_tform);

tampD = (max(d) < Tr);
geomD = (min(d_a) < Tr_a);

%geometric distortion detection

if (geomD)
    disp('Image_t has NOT undergone geometric manipulation.') ;
else
    disp('Image_t HAS undergone geometric manipulation.');
end

%tampering detection and localization

if (tampD)
    disp('Image_t is NOT tampered.') ;
else
    [blocks_RGB]=localization(d, d_a, blocks_RGB, Tr_a, Tr_1, Tr_2);
    disp('Image_t IS tampered.');
end

%% Show results

blocks_RGB=reshape(blocks_RGB,64,64);
blocks_RGB=cell2mat(blocks_RGB);
figure(1)
imshow(blocks_RGB)
authenticity = [geomD,tampD];
