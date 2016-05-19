function data = h5load(fname)
% Loads all variables in a hdf5 file into matlab

I = h5info(fname);
nd = length(I.Datasets);
for ii=1:nd
    eval(['data.' I.Datasets(ii).Name ' = ' ...
          'h5read(''' fname ''',''/' I.Datasets(ii).Name ...
          ''');']);
end

ng = length(I.Groups);
for ii=1:ng
    nd = length(I.Groups(ii).Datasets);
    gname = I.Groups(ii).Name;
    for jj=1:nd
    eval(['data.' I.Groups(ii).Datasets(jj).Name ' = ' ...
          'h5read(''' fname ''',''' gname '/' I.Groups(ii).Datasets(jj).Name ...
          ''');']);
    end
end

end