function [seqang,catang] = get_seqang(seqtlt,catmap,seqcat)
%  GET_SEQANG  Get sequence sample angles
%
%  Usage: [seqang,catang] = GET_SEQANG(seqtlt,catmap,seqcat)
%
%  where seqtlt - sequence sample tilts => in radians
%        catmap - angle/category mapping => 1 to 4
%        seqcat - sequence category => 1:orange or 2:blue
%
%  seqang contains the sample angles of the sequence (in rad), whereas catang
%  contains the categorization angle (in rad).
%
%  Valentin Wyart <valentin.wyart@ens.fr> - 09/2015

if nargin ~= 3
    error('Invalid list of input arguments!');
end

% create mapping matrix
matmap = zeros(8,2);
matmap(:,1) = (1:2:16)/16*pi; % orange w.r.t. angle/category mapping
matmap(:,2) = matmap(:,1)+pi/2; % blue w.r.t. angle/category mapping
matmap = mod(matmap,pi);

% apply mapping matrix
seqang = mod(seqtlt+matmap(catmap,seqcat),pi);
catang = mod(matmap(catmap,1)-pi/4,pi);

end