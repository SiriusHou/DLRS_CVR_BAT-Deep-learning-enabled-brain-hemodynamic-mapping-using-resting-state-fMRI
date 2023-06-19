function check = zach_contains(string,pattern)
    % takes a string and a pattern and checks to see if the pattern is
    % present in the string. If it is, check is returned as 1, otherwise 0.
    % (The function 'contains' is present in later versions of matlab; this
    % function was written for scripts being run in matlab 2013a)

    check = 0;
    for i = 1:length(string)-length(pattern)+1
        if strcmp(string(i:i+length(pattern)-1),pattern)
            check = 1;
            break
        end
    end

end