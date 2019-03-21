function s = spectrum(filePrefix, fileSuffix, fileNum)
    global datpath simpath;
    s = [];

    if isempty(fileNum)
        file = strcat(simpath, filePrefix,fileSuffix);
        s = load(file) * 6.2415091E9;
    else
        for i = fileNum
            file = strcat(datpath, filePrefix,num2str(i),fileSuffix);
            nrg = load(file) * 6.2415091E9;
            s = [s; nrg];
        end
    end
end
