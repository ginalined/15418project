clear all;
close all;
clc

% N = 100;
N = 10;
fov = 50;   % Define field of view of the camera in degrees

M = moviein(N);  % movie M
for number=1: N-1
    % extract data
    fname = sprintf('output/nbody_%04d.txt',number);
    obj_list = readfile(fname);
    
    % plot frame
    fig = figure();
    axis vis3d
    camva(fov);  % Set the camera field of view 
    campos([0 10 10]);
    camtarget([3, 3, 3]);
    camproj('perspective');
    camup([0, 1, 1]);

    for j = 1:size(obj_list, 1)
        % for each object, extract the triangle list
        triangle_list = obj_list.triangle_list;
        triangle_count = size(triangle_list, 2);
        C = rand(3);    % randomrize a color for the current object
        
        for i = 1:triangle_count
            % for each triangle, add the vertex to the plane vertex pool
            X = zeros(3);
            Y = zeros(3);
            Z = zeros(3);

            triangle = triangle_list(i);
            p1 = triangle.p1;
            p2 = triangle.p2;
            p3 = triangle.p3;

            X(1) = p1(1);
            X(2) = p2(1);
            X(3) = p3(1);

            Y(1) = p1(2);
            Y(2) = p2(2);
            Y(3) = p3(2);

            Z(1) = p1(3);
            Z(2) = p2(3);
            Z(3) = p3(3);

            % render triangle
            fill3(X, Y, Z, C);
            hold on;
        end
    end

    % collect frames
    M(number) = getframe(fig);
    size(M(number).cdata);
    close(fig);
end

output(M, 2);

close all;


% For each frame, extra a list of traingles with x,y,z positions
function obj_list = readfile(fname)
    obj_list = [];
    triangle_list = [];
    if exist(fname,'file')
        fid = fopen(fname);
        raw = fscanf(fid,'%f');
        fclose(fid);
    
        cur_obj = 0;
        for triangle_index = 1 : size(raw, 1)/10 - 1
            dp = raw(triangle_index * 10 + 1: (triangle_index + 1) * 10)
            triangle.p1 = [dp(1:3)];
            triangle.p2 = [dp(4:6)];
            triangle.p3 = [dp(7:9)];            
            
            if (dp(10) == cur_obj)
                triangle_list = [triangle_list, triangle];
            else
                obj.triangle_list = triangle_list;
                obj_list = [obj_list, obj];
                cur_obj = cur_obj + 1;
                triangle_list = [triangle];
            end
        end
        obj.triangle_list = triangle_list;
        obj_list = [obj_list, obj];   % wrap up
    else
        fprintf('File %s does not exist.\n', fname);    % throw exception
    end
end

% render or record the simulation video
function output(M, mode)
    if (mode == 1)
        % display the simulation movie
        movie(M, 10)
    else
        % record movie to avi
        v = VideoWriter('simulation.avi');
        v.FrameRate = 5;
        open(v);
        writeVideo(v,M);
        close(v);
    end
end
