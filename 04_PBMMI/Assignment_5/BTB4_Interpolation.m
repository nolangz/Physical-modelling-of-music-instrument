%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         Interpolation
%%%         Author:Nuolin Lai
%%%         Create Date:09/12/2020
%%%         Last modify date:26/03/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mtx = BTB4_Interpolation(p,a_q,fmode)
N = length(p);
%assert if N is an even number
if mod(N,2)~=0
    fprintf('please enter an even integer N.....');
    assert(mod(N,2)==0);
end
if fmode==1
    
    %calculate molecular
    aq  = repmat(a_q,1,N);
    aq_minor = aq-p;
    aq_prod = prod(aq_minor,2);
    aq_prod_rep = repmat(aq_prod,1,N);
    molecular = aq_prod_rep./aq_minor;

    %calculate demominator
    p_rep = repmat(p,N,1);
    p_rep_minor = -p_rep+p';
    p_eye = eye(N,N);
    denominator = prod(p_rep_minor+p_eye,2);

    %calculate mtx and fix error in special point
    mtx = molecular./denominator';
    mtx(1,N/2)=1;
    %do the same thing as above
elseif fmode==2

    % coefficient of the N Lagrange polynomials
    p   = -(N-1)/2:1:(N-1)/2;

    %calculate molecular
    aq  = repmat(a_q,1,N);
    aq_minor = aq-p;
    aq_prod = prod(aq_minor,2);
    aq_prod_rep = repmat(aq_prod,1,N);
    molecular = aq_prod_rep./aq_minor;

    %calculate demominator
    p_rep = repmat(p,N,1);
    p_rep_minor = -p_rep+p';
    p_eye = eye(N,N);
    denominator = prod(p_rep_minor+p_eye,2);

    %calculate mtx and fix error in special point
    mtx = molecular./denominator';
    mtx(isnan(mtx))=1;
    %plot the matrix
    plot(a_q,mtx,'LineWidth',1.5);
    legend(num2str((-(N-1)/2:(N-1)/2)'))
    xlabel('alpha');
    title('Lagrange Interpolation')
end
end