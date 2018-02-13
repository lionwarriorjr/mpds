function saveCsvAesFile(header, data, file, password )
%SAVEAESFILE Summary of this function goes here
%   Detailed explanation goes here
% save data to plainFile
plainFile = file;

[pathstr,name,ext] = fileparts(plainFile);
if ~isdir(pathstr)
    mkdir(pathstr);
end
csvwrite_with_headers(plainFile, data, header);


% encrypt plain file
if password ~= -1
    cmd = sprintf('aescrypt -e -p %s %s',password,plainFile) ;
    system(cmd)
    % delete plainFile
    delete(plainFile)
end



end