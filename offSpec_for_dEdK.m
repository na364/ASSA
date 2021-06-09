
function gammaOffSpec = offSpec_for_dEdK(Ei, dK)
% Returns the offSpec gamma equivalent to 'dK'

gammavec=22.2-55:0.025:22.2+55;
[dKvec,dkz] = k_xfer(beamprops('energy',Ei,3),gammavec,44.4);

for i=1:length(dK)
    indx=[]; threshold=0.01;
    while isempty(indx)
        indx=find(dKvec-dK(i)<threshold & dKvec-dK(i)>-threshold);
        threshold = threshold+0.01;
    end
    gammaOffSpec(i) = gammavec(indx(ceil(length(indx)/2)))-22.2;
end

end
