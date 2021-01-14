function options_struct = validate_classification_options_input(parsed_input, suppress_warnings)
%% validate input for classification parameters and output the options struct
if ~isfield(parsed_input,'opts_struct')
    options_struct = struct;
else
    options_struct = parsed_input.opts_struct;
end

%% get default setting for pairwise classification
if ~isfield(options_struct, 'pairwise') || isempty(options_struct)
    options_struct.pairwise = false; 
end

%% parameters specific to RSA and MCPA
if strcmp(func2str(parsed_input.test_handle), 'mcpa_classify') || strcmp(func2str(parsed_input.test_handle), 'rsa_classify')
    if ~isfield(options_struct, 'comparison_type')
        if isfield(options_struct, 'metric')
            switch options_struct.metric
                case {'spearman', 'pearson', 'kendall'}
                    if ~suppress_warnings
                        warning('Similarity or Dissimilarity space has not been defined. Based on your chosen metric, we will put the data in similarity space.')
                    end
                    options_struct.comparison_type = 'correlation';
                otherwise
                    if ~suppress_warnings
                        warning('Similarity or Dissimilarity space has not been defined. Based on your chosen metric, we will put the data in dissimilarity space.')
                    end                            
                    options_struct.comparison_type = 'distance';
            end
            
            switch options_struct.metric
                case '1-spearman'
                    options_struct.metric = 'spearman';
                case '1-pearson'
                    options_struct.metric = 'correlation';
            end

        else
            if ~suppress_warnings
                warning('Similarity or dissimilarity has not been defined. Putting the data in similarity space with the Spearman correlation.')
            end
            options_struct.comparison_type = 'correlation';
        end
    end
    if ~isfield(options_struct, 'metric')
        if strcmp(options_struct.comparison_type, 'correlation')
            options_struct.metric = 'spearman';
        else
            options_struct.metric = 'euclidean';
        end
    end
    
    if ~isfield(options_struct, 'tiebreak')
        options_struct.tiebreak = true;
    end
    
    if ~isfield(options_struct, 'exclusive')
        options_struct.exclusive = true;
    end
    
end

%% parameters for libsvm
if ~isfield(options_struct, 'libsvm_options') && strcmp(func2str(parsed_input.test_handle), 'libsvm_classify')
    options_struct.libsvm_options = '-s 0 -t 0 -q';
end
%% parameters specific to RSA
if ~isfield(options_struct, 'verbose') && strcmp(func2str(parsed_input.test_handle), 'rsa_classify')
    options_struct.verbose = 0;
end

end
