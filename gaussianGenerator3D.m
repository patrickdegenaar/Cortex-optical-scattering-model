% This function freely generates a gaussian matrix of set size and
% distribution. Typically the distribution will be 1/3 to 1/2 of the matrix
% size in order to make full use of the matrix box.

% I use an exponential decay function to determine the gaussian of the form: 
% exp(-(X^2+ y^2))

% Written by P. Degenaar on 19th March 2009

function gaussianImage = gaussianGenerator3D(size, distribution)

% create matrix
gaussianMatrix(1:size, 1:size) = 0;

% Scan through the matrix developing the elements
for i = 1:size
    for j = 1:size
        for k = 1:size
            x = i-ceil(size/2);
            y = j-ceil(size/2);
            z = k-ceil(size/2);
            
            gaussianMatrix(i,j,k) = 1*exp(-1 * (((x/distribution)^2)+ ((y/distribution)^2)+ ((z/distribution)^2)  )); 
        end
    end
end

% normalise
gaussianMatrix = gaussianMatrix./sum(sum(sum(gaussianMatrix)));

gaussianImage = gaussianMatrix;