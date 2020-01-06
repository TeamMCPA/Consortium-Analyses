
function [new_patterns, new_dimensions] = concatenate_dimensions(first_concat_dim, second_concat_dim, patterns)
%% This function creates a new folding dimension that is instances * subjects
% For when we want to do between subjects decoding but also want to
% maintain individual events
% Input: 
%   first_concat_dim = the first dimension we want in this concatenation.
%   second_concat_dim = the second dimension we want in this concatenation.
%   patterns = MCPA pattern matrix in the form time x condition x channel x
%   repetition x session x subject

% Output: 
%   new_patterns = the new pattern matrix with the concatenated
%   dimensions. The new dimension will be in the place of the
%   first_concat_dim and will be the size of
%   first_concat_dim*second_concat_dim
%   new_dimensions = a cell array of the new dimensions. The concatenated
%   dimension is represented as 'first_concat_dim+second_concat_dim'.

try
    if first_concat_dim == 2 && second_concat_dim == 4
        %% combine dimensions 2 and 4
        new_patterns = nan(size(patterns,1),...
                        (size(patterns,2)*size(patterns,4)),...
                        size(patterns,3),...
                        size(patterns,5),...
                        size(patterns,6));

        to_add = size(patterns,2);

        for dim = 1:size(patterns,4)
            start_idx = (to_add*dim) - to_add + 1;
            stop_idx = to_add * dim;
            new_patterns(:,start_idx:stop_idx,:,:,:) = squeeze(patterns(:,:,:,dim,:,:)); 
        end
        new_dimensions = {'time', 'condition+repetition','channel', 'session', 'subject'};

    elseif first_concat_dim == 2 && second_concat_dim == 5
        %% combine dimensions 2 and 5
        new_patterns = nan(size(patterns,1),...
                        (size(patterns,2)*size(patterns,5)),...
                        size(patterns,3),...
                        size(patterns,4),...
                        size(patterns,6));

        to_add = size(patterns,2);

        for dim = 1:size(patterns,5)
            start_idx = (to_add*dim) - to_add + 1;
            stop_idx = to_add * dim;
            new_patterns(:,start_idx:stop_idx,:,:,:) = squeeze(patterns(:,:,:,:,dim,:)); 
        end
        new_dimensions = {'time', 'condition+session', 'channel', 'repetition', 'subject'};

    elseif first_concat_dim == 4 && second_concat_dim == 5
        %% combine dimensions 4 and 5
        new_patterns = nan(size(patterns,1),...
                        size(patterns,2),...
                        size(patterns,3),...
                        (size(patterns,4)*size(patterns,5)),...
                        size(patterns,6));

        to_add = size(patterns,4);

        for dim = 1:size(patterns,5)
            start_idx = (to_add*dim) - to_add + 1;
            stop_idx = to_add * dim;
            new_patterns(:,:,:,start_idx:stop_idx,:) = squeeze(patterns(:,:,:,:,dim,:)); 
        end
        new_dimensions = {'time', 'condition', 'channel', 'repetition+session', 'subject'};
    elseif first_concat_dim == 4 && second_concat_dim == 6
        %% combine dimensions 4 and 6
        new_patterns = nan(size(patterns,1),...
                        size(patterns,2),...
                        size(patterns,3),...
                        (size(patterns,4)*size(patterns,6)),...
                        size(patterns,5));

        to_add = size(patterns,4);

        for dim = 1:size(patterns,6)
            start_idx = (to_add*dim) - to_add + 1;
            stop_idx = to_add * dim;
            new_patterns(:,:,:,start_idx:stop_idx,:) = squeeze(patterns(:,:,:,:,:,dim)); 
        end
        new_dimensions = {'time', 'condition', 'channel', 'repetition+subject', 'session'};
    elseif first_concat_dim == 5 && second_concat_dim == 6
        %% combine dimensions 5 and 6 
        new_patterns = nan(size(patterns,1),...
                        size(patterns,2),...
                        size(patterns,3),...
                        size(patterns,4),...
                        (size(patterns,5)*size(patterns,6)));

        to_add = size(patterns,5);            

        for dim = 1:size(patterns,6)
            start_idx = (to_add*dim) - to_add + 1;
            stop_idx = to_add * dim;
            new_patterns(:,:,:,:,start_idx:stop_idx) = squeeze(patterns(:,:,:,:,:,dim)); 
        end
        new_dimensions = {'time', 'condition', 'channel', 'repetition', 'session+subject'};
    else
        fprintf('Those dimensions cannot be concatenated.');
    end
catch
    error('Cannot support this form of concatenation. Be sure pattern dimensions are in form of time x condition x channel x repetition x session x subject.')
end
    
end













