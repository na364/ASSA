
function dK = dK_for_gammaOfSpec_Ei(Ei, gammaOffSpec, tot_angle)

% Returns the dK for offSpec gamma
if ~exist('tot_angle','var'), tot_angle=44.4;end
[dK,dkz] = k_xfer(beamprops('energy',Ei,3),gammaOffSpec+tot_angle/2,tot_angle);


