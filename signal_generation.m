%{
function for signal generation
T -> N x L
Phi -> N x M
W -> M x L
En -> N x L (noise)
N << M
%}

function T = signal_generation(Phi, W, En)
    T = Phi*W + En;
end


    
