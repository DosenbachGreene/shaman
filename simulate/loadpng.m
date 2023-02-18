function img = loadpng(file)
% Load a portable network graphic file, scale to 1, and mirror upper
% triangular portion about the diagonal.

% Read in the image file.
img = double(imread(file));

% Convert to monochrome.
img = mean(img,3);

% Scale image data to range between 0 and 1.
img = (double(img) - 128) ./ 128;

% Delete the lower triangular part of the image.
img = img - tril(img);

% Reflect image of upper triangular part about diagnoal to create
% correlation matrix.
img = eye(size(img)) + 0.5 .* (img + img');

end