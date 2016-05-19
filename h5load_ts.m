function data = h5load_ts(fname,ti)
% Loads all variables in a hdf5 file into matlab

    data.ti = ti;
I = h5info(fname);
nd = length(I.Datasets);
for ii=1:nd
        sz = I.Datasets(ii).Dataspace.Size;
        if (length(sz) == 4)
            eval(['data.' I.Datasets(ii).Name ' = ' ...
                  'h5read(''' fname ''',''/' I.Datasets(ii).Name ...
                  ''',[1 1 1 ti],[sz(1) sz(2) sz(3) 1]);']);        
        else
            eval(['data.' I.Datasets(ii).Name ' = ' ...
                  'h5read(''' fname ''',''/' I.Datasets(ii).Name ...
                  ''');']);
        end
end

ng = length(I.Groups);
for ii=1:ng
    nd = length(I.Groups(ii).Datasets);
    gname = I.Groups(ii).Name;
    for jj=1:nd
        sz = I.Groups(ii).Datasets(jj).Dataspace.Size;
        if (length(sz) == 4)
            eval(['data.' I.Groups(ii).Datasets(jj).Name ' = ' ...
                  'h5read(''' fname ''',''' gname '/' I.Groups(ii).Datasets(jj).Name ...
                  ''',[1 1 1 ti],[sz(1) sz(2) sz(3) 1]);']);        
        else    
            eval(['data.' I.Groups(ii).Datasets(jj).Name ' = ' ...
                  'h5read(''' fname ''',''' gname '/' I.Groups(ii).Datasets(jj).Name ...
                  ''');']);
        end
    end
end

end