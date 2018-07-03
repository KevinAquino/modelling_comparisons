% Function here to generate the time series using the time series from
% averaging inside the parcel -- specifically for HCP nothing else just yet
% -- I might make it more sophisticated later for other pracellations

function out = parcelTimeSeries(timeSeries,parcellation)
    [~, label, ctab] = read_annotation(parcellation);
    
    % Initializing the matrix.
    out = zeros(length(ctab.struct_names)-1,size(timeSeries,2));
    
%     Looping through the array to take the average of each parcel
    for parcel=2:length(ctab.struct_names)
        inds = find(label==ctab.table(parcel,5));
        sampledSeries = timeSeries(inds,:);
        out(parcel-1,:) = mean(sampledSeries,1);
    end

end