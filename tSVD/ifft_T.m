function [U,S,V] = ifft_T(U,S,V,la)
for i = la:-1:3
     U= ifft(U,[],i);
     S= ifft(S,[],i);
     V= ifft(V,[],i);
end
end