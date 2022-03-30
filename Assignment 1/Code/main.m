clc; clear;
close all;
addpath("./MMread/");
T = 3;
mean_noise = 0;
sigma = 2;
p = 8;
% Main frames
path = '../cars.avi';
% Reading video from the folder, coverting to grayscale and extrating T = 3 frames
video = mmread (path,1:T);
frames_array = zeros(video.height, video.width, T);
for k = 1 : T
    frames_array(:,:,k) = im2double(rgb2gray(video.frames(k).cdata));
end

%changing the frame to grayscale frame
frames_array = frames_array(end-120:end, end-240:end, :);


figure();
montage([frames_array(:,:,1), frames_array(:,:,2), frames_array(:,:,3)]);
title("Original frames");

f_reconst_3 = video_process(frames_array, T ,p,mean_noise, sigma);
rmse(f_reconst_3, frames_array)

T = 5
% Reading video from the folder, coverting to grayscale and extrating T = 3 frames
video = mmread (path,1:T);
frames_array2= zeros(video.height, video.width, T);
for k = 1 : T
    frames_array2(:,:,k) = im2double(rgb2gray(video.frames(k).cdata));
end

%changing the frame to grayscale frame
frames_array2 = frames_array2(end-120:end, end-240:end, :);

figure();
montage([frames_array(:,:,1), frames_array(:,:,2), frames_array(:,:,3)]);
title("Original frames");

f_reconst_5 = video_process(frames_array2, 5 ,p,mean_noise, sigma);
rmse(f_reconst_5, frames_array2);

T = 7
% Reading video from the folder, coverting to grayscale and extrating T = 3 frames
video = mmread (path,1:T);
frames_array3 = zeros(video.height, video.width, T);
for k = 1 : T
    frames_array3(:,:,k) = im2double(rgb2gray(video.frames(k).cdata));
end

%changing the frame to grayscale frame
frames_array3 = frames_array3(end-120:end, end-240:end, :);


figure();
montage([frames_array3(:,:,1), frames_array3(:,:,2), frames_array3(:,:,3)]);
title("Original frames");

f_reconst_7 = video_process(frames_array3, T ,p,mean_noise, sigma);
rmse(f_reconst_7, frames_array3);

T = 5
path = '../flame.avi'
% Reading video from the folder, coverting to grayscale and extrating T = 3 frames
video = mmread (path,1:T);
frames_array4 = zeros(video.height, video.width, T);
for k = 1 : T
    frames_array4(:,:,k) = im2double(rgb2gray(video.frames(k).cdata));
end

%changing the frame to grayscale frame
frames_array4 = frames_array4(end-120:end, end-240:end, :);

figure();
montage([frames_array4(:,:,1), frames_array4(:,:,2), frames_array4(:,:,3)]);
title("Original frames");

flame_reconst = video_process(frames_array4, T ,p,mean_noise, sigma);
rmse(flame_reconst, frames_array4);

function f_reconst = video_process(frames_array, T, p, mean_noise , sigma) 
    height = size(frames_array, 1);
    width = size(frames_array, 2);

    S  = randi([0,1], height, width, T);
    codedsnapshot = S .* frames_array;
    I = sum(codedsnapshot,3);
    
    % Adding noise to the codedsnapshot
    noise = normrnd(mean_noise, sigma/255, height, width);
    I_final = I + noise;
    figure();
    imshow(mat2gray(I_final));
    title("Codedsnapshot");
    axis tight;

    dct = dctmtx(p);
    psi = kron(dct, dct);

    %Taking same psi for three frames.
    psi_t = kron(eye(T), psi);
    f_reconst = zeros(height, width, T);
    count = zeros(height, width, T);

    %patch wise reconstruction
    for i = 1:height - p + 1
        for j = 1: width - p + 1
            
            % \Phi_t is a diagonal matrix and f_t is 
            % the t-th vectorized frame.
            % y has n pixels, and so does each 
            % f_t and each \Phi_t is a n x n diagonal matrix. 
            patch = I(i:i+p-1, j:j+p-1);
            patch_S = S(i:i+p-1, j:j+p-1,:);

            % \Phi_{it} is a p^2 times p^2 diagonal matrix
            phi = [];
            for t = 1:T
                st = patch_S(:, :, t);
                phi = [phi diag(reshape(st, p^2, 1))];
            end
            
             A = phi * psi_t';
            
    
            % f_{it} is a p x p patch from the t-th frame 
            % but expressed as a vector of size p^2 by 1
            % y_i is also a p^2 times 1 vector.
            y = reshape(patch, p^2, 1);
            
            % theta_i is a vector of size Tp^2 by 1 
            % obtained by concatenating all T vectors theta_{it}
            [theta, ~]= OMP(A, y,.035 , 128);
            I_t = psi_t' * theta;
            for t = 1: T
                f_t = reshape(I_t((t-1)*p^2+1:t*p^2), p,p);
%                 f_t = idct2(f_t);
%                 f_t = f_t*;
                f_reconst(i:i+p-1, j: j+p-1, t) = f_reconst(i:i+p-1, j: j+p-1, t) + f_t;

                % Repeat the reconstruction for all overlapping 
                % patches and average across the overlapping 
                % pixels to yield the final reconstruction.
                count(i:i+p-1, j: j+p-1, t) = count(i:i+p-1, j: j+p-1, t) + ones(p, p);
            end
            
        end
    end
    
    % averaging across the overalapping stages . 
    % Using of elementwise counts to take average for each pixels.
    f_reconst = f_reconst ./ count;

    % Printing frames
    for  t = 1:T
         figure();
         montage([mat2gray(f_reconst(:,:,t)), frames_array(:,:,t)]);
        title("T = "+ t);
    end

end

% Subtracting values from each frame for getting all error values and
% taking the mean.
function rmse = rmse(f_reconst, E)
    rmse = mean((f_reconst - E).^2, 'all')/mean(E.^2, 'all');
    fprintf("rmse = %.3f \n", rmse);
end

