
function parsed_opts = parse_opts(opts)
    if isempty(opts)
        opts = struct;
    end
    % get the fields and values from the opts struct
    fields = fieldnames(opts)';
    values = struct2cell(opts)';
    
    input = cell(1,2*length(fields)); % initialize an empty cell array
    
    
    % combine the struct fields followed by their respective values
    for i = 1:length(input)
        if rem(i,2) == 0
            input(i) = values(i/2); % set even indices to the struct values
        else
            input(i) = fields(ceil(i/2)); % set odd indices to the struct field names
        end
    end

    parsed_opts = input;
    
end