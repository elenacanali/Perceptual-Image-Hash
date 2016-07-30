function [F_0, Gs1, Gs2]=extract_features(image_0, s1, s2, tc,f_tform)
% Input: 
% image_0 (of class single) : the image to hash
% s1, s2 : the projection rates for the compression matrices
% tc: a String indicating the type of compression matrix to use
% -> 'Bernoulli' or 'Gauss'

%% Parameters initialization

% initialize intermediate hash F_0:={GT, GB, S(x,y), M, N}
F_0=cell(1,5);
[F_0{4}, F_0{5}]=size(image_0);

% initialize the dimension of the features
features_dimension = 0;

%% (a) KEY-POINT-BASED FEATURES 
% Extract the SIFT features or the SURF features of the image and apply a
% 1-level db1 wavelet transform to each them separately.

if strcmp(f_tform, 'SURF')
    
    % Detect the SURF features of the image and extract them, alongside
    % with their location, then cast the features from single to double and
    % saves the coordinates, deleting all additional information.
    
    points = detectSURFFeatures(image_0);
    [features,locations] = extractFeatures(image_0,points);
    F_0{1}=double(features'); % cast features for future computation
    F_0{3}=double(locations.Location);

    for i=1:size(F_0{1},2)
       % apply 1-level db1 transformation to F{1}
       [F_0{1}(1:32, i), F_0{1}(33:64, i)]=dwt(F_0{1}(:,i), 'db1');
    end
    
    % set dimension of the features
    features_dimension = 64;
end
if strcmp(f_tform, 'SIFT')
                  
    % Compute the SIFT frames (keypoints) and descriptors. The matrix S has
    % a column for each frame. A frame is a disk of center S(1:2), scale S(3)
    % and orientation S(4). Since we need only S(1:2), delete other components.

    [F_0{3}, F_0{1}] = vl_sift(image_0);
    F_0{1}=double(F_0{1}); % cast F{1} for future computation
    F_0{3}=F_0{3}(1:2,:)';

    for i=1:size(F_0{1},2)
        % apply 1-level db1 transformation to F{1}
        [F_0{1}(1:64, i), F_0{1}(65:128, i)]=dwt(F_0{1}(:,i), 'db1');
    end 
    
    % set dimension of the features
    features_dimension = 128;
end
                  
%% (b) BLOCK-BASED FEATURES
% Divide a single precision image into 8x8 non overlapping
% blocks and save it in to an array of n_blocks=(MxN/P^2) blocks

n_blocks=floor(F_0{4}/8)*floor(F_0{5}/8);
blocks=mat2cell(image_0, 8*ones(1,F_0{4}/8), 8*ones(1,F_0{5}/8));
blocks=reshape(blocks,1,n_blocks);
                  
% For Each block, compute DCT coefficients
CB=zeros(n_blocks, 8, 8);
for i=1:n_blocks
   CB(i, :, :)=dct2(blocks{i}); 
end
load watson; % t=Watson's matrix 
F_0{2}=zeros(64, n_blocks); % Weighted dct matrix (block wise) and zigzag

for i=1:n_blocks   
    tmp=reshape(CB(i,:,:), 8, 8); % need by ./t
    tmp=tmp./t;
    F_0{2}(:, i)=zigzag(tmp)';
end

%% (c) COMPRESSION AND PROJECTION

% create Gaussian random matrices Gs1, Gs2
if strcmp(tc,'Gauss')
    Gs1=randn(ceil(features_dimension*s1), features_dimension);
    Gs2=randn(ceil(64*s2), 64);
end
% alternative: create bernoullian matrices
if strcmp(tc,'Bernoulli')
    Gs1 = randi([0,1],ceil(features_dimension*s1),features_dimension) * 2 - 1;
    Gs2 = randi([0,1],ceil(64*s2),64) * 2 - 1;
end

% Apply compression
F_0{1}=Gs1*F_0{1};
F_0{2}=Gs2*F_0{2};

