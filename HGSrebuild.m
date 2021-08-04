function [species,n0] = HGSrebuild(species,n0)
%**************************************************************************
%
% [species,n0] = HGSrebuild(species,n0)
%
%**************************************************************************
%
% HGSrebuild transform the mixtures in to their species. It also transforms
% the mols of the mixtures with the proportions assigned.
%
%**************************************************************************
% Inputs:
%--------------------------------------------------------------------------
% species --> String or numbers of species
% n0 --> [mols] Mols
%
% Outputs:
%--------------------------------------------------------------------------
% species --> String of species
% n0 --> [mols] Mols
%
%**************************************************************************
% *HGS 2.0
% *By Caleb Fuster, Manel Soria and Arnau Miró
% *ESEIAAT UPC    



global HGSdata; HGSload;

if isnumeric(species)
    for ii=1:length(species)
        if species(ii)<=length(HGSdata.Mm)
            nspec{ii} = HGSdata.name{species(ii)};
        else
            nspec{ii} = HGSdata.comb{species(ii)-length(HGSdata.Mm)};
        end
    end
    species = nspec;
end

if ischar(species)
   species = {species};
end


jj = 1;
while jj <= length(species)
   id = HGSid(species{jj});
   if id > length(HGSdata.Mm)
       for ii=1:length(HGSdata.comb)
           if strcmp(HGSdata.comb{ii},species{jj})
                newspecies = species(1,1:jj-1);
                newn = n0(1,1:jj-1);              

                for kk=1:length(HGSdata.cspec{ii})
                    newspecies(jj+(kk-1)) = HGSdata.cspec{ii}(kk);  
                    newn(jj+(kk-1)) = n0(jj)/100*HGSdata.cper{ii}(kk);
                end

                species = [newspecies, species{1,jj+1:end}];
                n0 = [newn, n0(1,jj+1:end)];

           end
       end
   else
       jj = jj+1;
   end   
end

Elem = unique(species);
if length(Elem) ~= length(species)
    n = zeros(1,length(Elem));
    for ii = 1:length(species)
        for jj = 1:length(Elem)
            if strcmp(Elem{jj},species{ii})
                n(jj) = n(jj) + n0(ii);
            end
        end  
    end
    n0 = n;
    species = Elem;
end

end