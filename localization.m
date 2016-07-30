function [blocks_RGB]=localization(d, d_a, blocks_RGB, Tr_a, Tr_1, Tr_2)
%% Step (6) and Step (7)
% Estimate the Euclidean distance d_a between each pair of matched
% feature points in set S_0 and S_t. The similarity of each block
% is compared to Tr_1 (if no geometric manipulation is present) or
% to Tr_2 (otherwise) and the block is marked if tampered

n_blocks=size(blocks_RGB,2);
if min(d_a)>Tr_a % testing for geometric transformation
    for i=1:n_blocks
        if d(i)>Tr_2
            blocks_RGB{i}(:,:,1)= ones(8); % if compromised, red
            blocks_RGB{i}(:,:,2)=zeros(8);
            blocks_RGB{i}(:,:,3)=zeros(8);
        end
    end
else
    for i=1:n_blocks % generic compromission, Tr_2 > Tr_1
        if d(i)>Tr_1
            blocks_RGB{i}(:,:,1)= ones(8); % if compromised, red
            blocks_RGB{i}(:,:,2)=zeros(8);
            blocks_RGB{i}(:,:,3)=zeros(8);
        end
    end
end