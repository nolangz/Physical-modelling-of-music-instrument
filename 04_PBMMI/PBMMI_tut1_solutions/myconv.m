% -------------------------------------------------
% FIR CONVOLUTION EQUATION
% Alistair Carson 2021
% -------------------------------------------------
function y = myconv(h,x)

% switch h and x for speed if necessary
if length(h) > length(x)
    temp = x;
    x = h;
    h = temp;
end

% lengths
Nx = length(x);
Nb = length(h)-1;

% initialise
y = zeros(Nx + Nb,1);

% convolution
% for m = 0:Nb
%     y = y + h(m+1)*[zeros(m,1); x; zeros(Nb-m, 1)];
% end

% alternative, faster method
xpad = [zeros(Nb, 1); x; zeros(Nb,1)];
for m = 0:Nb
    y = y + h(m+1)*xpad(Nb+1-m:end-m);
end


end
