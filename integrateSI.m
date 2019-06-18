function [ int_si ] = integrateSI( bksub_si, en, intwin )
%integrateSI Integrates signal in background-subtracted spectrum image
% inputs:
%   bksub_sisi -- the background-subtracted spectrum image, ordered [x,y,en]
%   en -- a vector of energy values corresponding to the energy axis.
%   intwin -- energy values specifying energy range for integration, in
%             format [enmin, enmax]
% outputs:
%   int_si -- 2D image of integrated signal

log_en = en>=min(intwin) & en<=max(intwin);
p = size(bksub_si);

if length(p)==3 % you have 3d si
    int_si = sum(bksub_si(:,:,log_en),3);
else
    if length(en) == p(2) % energy is 2nd dim
        int_si = sum(bksub_si(:,log_en),2);
    elseif length(en) == p(1) % energy is 1st dim
        int_si = sum(bksub_si(log_en,:),1);
    end
end

end

