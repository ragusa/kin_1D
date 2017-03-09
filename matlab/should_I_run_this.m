function true_false = should_I_run_this(list,name)


true_false = false;

for i=1:length(list)
    true_false = strcmpi(list{i},name);
    if true_false
        break
    end
end

return
end

