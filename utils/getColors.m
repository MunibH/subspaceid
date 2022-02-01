function clrs = getColors()

reds = {[128 0 22]
        [160 0 28]
        [192 0 33]
        [255 0 43]};
    
reds = cellfun(@(x) x./255, reds,'UniformOutput',false);
    
blues = {[0 4 58]
        [0 41 98]
        [0 78 137]
        [64 123 167]};

blues = cellfun(@(x) x./255, blues,'UniformOutput',false);

clrs.reds = reds;
clrs.blues = blues;

end