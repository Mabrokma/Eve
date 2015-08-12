%Make an expression prediction based on a set of model parameters and
%transcription factor concentrations
function expr = useThermoModel(Params, sites, tfConc)
traceLength = size(tfConc, 1);
nBins = size(tfConc,2);
nSites = length(sites);

%loop over all time and space
for t = 1:traceLength
    for x = 1:nBins
        fracOcc = NaN(1,length(sites));
        %calculate fractional occupancy of each site
        for i = 1:nSites
            m_i = sites(i).m; %Start of current site
            n_i = sites(i).n; %end of current site
            conc = tfConc(sites(i).Name); %profile of tf
            
            top = Params.E(sites(i).Name)*conc(t,x);
            bottom = 1 + top;
            
            %Extract overlapping sites
            overlappingSitesIndices = ...
                (([sites.n] > m_i) & ([sites.n] < n_i))...
                | (([sites.m] > m_i) & ([sites.m] < n_i));
            overlappingSites = sites(overlappingSitesIndices);
            
            for j = 1:length(overlappingSites)
                conc = tfConc(ovelappingSites(j).Name);
                bottom = bottom ...
                    + Params.E(overlappingSites(j).Name)*conc(t,x);
            end
            
            fracOcc(i) = top / bottom;
        end
        
        %correct for quenching
        
        %Calculate number of recruited adapter factors
        
        
        %Calculate rate of mRNA production
        
        %THIS IS WHERE OUR DATA IS DIFFERENT FROM THEIRS.  
        %INSTEAD MAYBE CALCULATE PROBABILITY OF BEING ON?
    end
end