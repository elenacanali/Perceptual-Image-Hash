function [d, d_a, blocks_RGB]=verifier(image_t, F_E, K, Gs1, Gs2, Tr_a,f_tform)
%% Step (1) 
% Decrypting
F_D=decrypt(F_E, K);
image_t = single(image_t);

%% Step (2)

%  Key point-based features extraction from image_t

F_t=cell(1,5);
F_t{4}=F_D{4};
F_t{5}=F_D{5};

if strcmp(f_tform, 'SURF')
    
    % Detect the SURF features of the image and extract them, alongside
    % with their location, then cast the features from single to double and
    % saves the coordinates, deleting all additional information.
    
    points = detectSURFFeatures(image_t);
    [features,locations] = extractFeatures(image_t,points);
    F_t{1}=double(features'); % cast features for future computation
    F_t{3}=double(locations.Location);

    for i=1:size(F_t{1},2)
       % apply 1-level db1 transformation to F{1}
       [F_t{1}(1:32, i), F_t{1}(33:64, i)]=dwt(F_t{1}(:,i), 'db1');
    end
end

if strcmp(f_tform, 'SIFT')
                  
    % Compute the SIFT frames (keypoints) and descriptors. The matrix S has
    % a column for each frame. A frame is a disk of center S(1:2), scale S(3)
    % and orientation S(4). Since we need only S(1:2), delete other components.
    
    [F_t{3}, F_t{1}] = vl_sift(image_t);
    F_t{1}=double(F_t{1}); % cast F{1} for future compuation
    F_t{3}=F_t{3}(1:2,:)';

    for i=1:size(F_t{1},2)
        % apply 1-level db1 transformation to F{1}
        [F_t{1}(1:64, i), F_t{1}(65:128, i)]=dwt(F_t{1}(:,i), 'db1');
    end 
end
 
F_t{1}=Gs1*F_t{1}; % apply Gs1 projection

%% Step (3)
% Use the proper matching algorithm to obtain n pairs of most similar
% features and relative coordinates, with 4 <= n <= min(m, m_t, 40)

if strcmp(f_tform,'SURF')
    matches=matchFeatures(F_D{1}',F_t{1}','Unique',true);
end
if strcmp(f_tform,'SIFT')
    matches=vl_ubcmatch(F_D{1}, F_t{1})';
end
                  
% Update the feature point sets as matched feature points
S_0 = F_D{3}(matches(:,1),:);
S_t = F_t{3}(matches(:,2),:);

%% Step (4 and 7)
% Detect and estimate geometric manipulation.
% Remark. Step (7) is now executed here to avoid distortion of
% authentic images due to numerical imprecisions
                  
% Estimate the Euclidean distance d_a between each pair of matched
% feature points in set S_0 and S_t.

d_a=zeros(1,size(S_0,1));
for i=1:size(d_a,2)
    d_a(i)=norm(S_0(i,:)-S_t(i,:)); 
end

% If geometric mainpulation is detected, estimate its associated
% affine transform PI and rectify the image
                  
if(min(d_a)>Tr_a)
                    
    % Define the linear system given by the matching point
    M_0=ones(3,size(S_0,1));
    M_t=M_0;
    M_0(1,:)=S_0(:, 1)'; % y coordinates of set S_0
    M_0(2,:)=S_0(:, 2)'; % x coordinates of set S_0
    M_t(1,:)=S_t(:, 1)'; % y coordinates of set S_t
    M_t(2, :)=S_t(:, 2)'; % x coordinates of set S_t
    
    % Solve the linear system using EM algorithm.
    % This will help to exclude strong outliers, i.e. false positives
    % in the matching algorithm.
    
    % Initial conditions
    sigma = 10000; % variance on Gaussian form, used big number to avoid 
                   % nullification of weights with early values.
    PI = diag(ones(1,3)); 

    count = 0;
    R = zeros(1,size(M_0,2));
    WM0 = zeros(3,size(M_0,2));
    WMT = zeros(3,size(M_0,2));
    while(1) 

        % E-STEP
        % calculate residual
        res = (PI*M_0 - M_t);
        res = res(1:2,:); % the third line is forced to be all ones.
        for i = 1 : size(res,2);
            R(i) = norm(res(:,i));  
        end
        
        % compute weight for each point.
        E = exp(-R.^2/sigma); 

        % M-STEP

        % Use weighted least squares to solve the linear system
        for i = 1 : size(M_0,2)

            WM0(:,i) = E(i)*M_0(:,i);
            WMT(:,i) = E(i)*M_t(:,i);

        end

        % Rem. linsolve uses least squares to solve overdetermined linear systems
        PInew = linsolve(WM0',WMT');
        PInew(:,3) = [0,0,1]'; % force PI to be affine relation;

        if(sum(abs(PI(1:2,:))) == 0)
            PI(1,1) = 1;
            PI(2,2) = 1;
        end

        d = 0;
        for i = 1 : 2
            for j = 1 : 3
                d = d + abs(PI(i,j)-PInew(i,j));
            end
        end

        if (d < 0.001) % stop condition
            break;
        end

        PI = PInew;
        count = count+1;

         if (count > 20) % escape for some non-convergent cases
             break;
         end
 
         sigma = sum(E.*R.^2)/sum(E); % update sigma
     end
    
    PI = inv(PI);
    PI(:,3) = [0,0,1]'; % remove errors due to numerical precision,
                        % forcing PI to be affine
    if(max(PI(1,1),PI(2,2)) > 5 )
        warning('Risky affine transform miscalculation');
        PI = diag(ones(3,1));
    end
    tform = affine2d(PI);

    % apply the affine transform and trim the edges left by the transformation

    image_t = imwarp(image_t,tform);
    [h,w]= size(image_t);
    originalH = 8*round(F_t{4}/8);
    originalW = 8*round(F_t{5}/8);
    image_t = imcrop(image_t,[floor(0.5*(w-originalW)),floor(0.5*(h-originalH)),originalW,originalH]);

end

image_t = imresize(image_t,[8*round(F_t{4}/8), 8*round(F_t{5}/8)]);
%% Step (5)
% Block-based features extraction

n_blocks = round(F_t{4}/8)*round(F_t{5}/8);
dimension=8*ones(1,64);
blocks_t=mat2cell(image_t, dimension, dimension);
blocks_RGB=mat2cell(image_t(:,:,[1 1 1]),dimension,dimension,3);
blocks_t=reshape(blocks_t,1, n_blocks);
blocks_RGB=reshape(blocks_RGB,1, n_blocks);
CB_t=zeros(n_blocks, 8, 8);

for i=1:n_blocks
   CB_t(i, :, :)=dct2(blocks_t{i}); 
end

load watson;

F_t{2}=zeros(64, n_blocks);
for i=1:n_blocks
    tmp=reshape(CB_t(i,:,:), 8, 8);
    tmp=tmp ./ t;
    F_t{2}(:, i)=zigzag(tmp)';
end
F_t{2}=Gs2*F_t{2}; % apply Gs2 projection



%% Step (6)
% measure similarity between image_0-blocks and image_t-blocks
% using Euclidean distance

d=zeros(1, n_blocks);
for i=1:n_blocks
    d(i)=norm(F_D{2}(:,i)-F_t{2}(:,i));
end

