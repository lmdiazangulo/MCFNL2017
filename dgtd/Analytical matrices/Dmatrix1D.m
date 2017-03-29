function [ D ] = DMatrix1D(N,sDir)
    D = Pmatrix1D(N,sDir) * DmatA(N,1) * Pmatrix1D(N,sDir);
end

    
