function options_struct = validate_classification_options_input(MCP_struct,parsed_input, suppress_warnings)
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


%% MCPA only parameters
if strcmp(func2str(parsed_input.test_handle), 'mcpa_classify') 
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

%% RSA only parameters
if strcmp(func2str(parsed_input.test_handle), 'rsa_classify')
    if ~isfield(options_struct, 'weighted_average')
        options_struct.weighted_average = false;
    end
    
    if options_struct.weighted_average
        for subj = 1:length(MCP_struct)
            for sess = 1:length(MCP_struct(subj).Experiment.Runs)
                session_inds = MCP_struct(subj).Experiment.Runs(sess).Index;
                trials_per_session(subj, sess) = length(find(MCP_struct(subj).fNIRS_Data.Onsets_Matrix(session_inds,:)));  
            end
        end         
        trials_per_session(trials_per_session==0) = NaN;
        options_struct.trials_per_session = trials_per_session;    
    end
    
    if ~isfield(options_struct, 'similarity_function')
        if isfield(options_struct, 'metric')
            switch options_struct.metric
                case {'spearman', 'pearson', 'kendall'}
                    if ~suppress_warnings
                        warning('Similarity or Dissimilarity space has not been defined. Based on your chosen metric, we will put the data in similarity space.')
                    end
                    options_struct.similarity_function = @create_similarity_matrix;
                otherwise
                    if ~suppress_warnings
                        warning('Similarity or Dissimilarity space has not been defined. Based on your chosen metric, we will put the data in dissimilarity space.')
                    end                            
                    options_struct.similarity_function = @create_dissimilarity_matrix;
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
            options_struct.similarity_function = @create_similarity_matrix;
        end
    end
    
     if ~isfield(options_struct, 'metric')
        if strcmp(func2str(options_struct.similarity_function), 'create_similarity_matrix')
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
