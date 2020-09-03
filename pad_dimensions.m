function padded_inds = pad_dimensions(current_dimensions, dim_to_match, value_to_insert)
%% create a cell array that will flexibly subset out one dimension out of a multidimensional matrix
% current_dimensions: list of the dimensions your matrix currently has
% (e.g. {'time', 'condition', 'feature'})
% dim_to_match: dimension you want to subset out (e.g. 'time')
% value_to_inset: what indices of that dimension do you want to subset
% (e.g. 1:5)

padded_inds = repmat({':'},1,length(current_dimensions)); % Create index structure with all-elements in all-dimensions
padded_inds{strcmp(current_dimensions,dim_to_match)} = value_to_insert;% In whichever dimension matches 'session', substitute the incl_sessions vector

end



