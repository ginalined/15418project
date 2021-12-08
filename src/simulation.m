clear all;
close all;
clc

% Nframe = 100;
Nframe = 10;
fov = 50;   % Define field of view of the camera in degrees

M = moviein(Nframe);  % movie M
for number=1: Nframe
    % plot frame
    fig = figure();
    % axis vis3d
    % camva(fov);  % Set the camera field of view 
    % campos([0 10 10]);
    % camtarget([3, 3, 3]);
    % camproj('perspective');
    % camup([0, 1, 1]);

    fname = sprintf("nbody_%04d.inp", number)
    obj_list = readTriangle(fname);
    perFrameRender(obj_list);

    % collect frames
    M(number) = getframe(fig);
    size(M(number).cdata);
    close(fig);
end

output(M, 2);

close all;


% For each frame, extra a list of traingles with x,y,z positions
function obj_list = readTriangle(fname)
    obj_list = [];
    triangle_list = [];
    
    if exist(fname,'file')
        fid = fopen(fname);
        raw = fscanf(fid,'%f');
        TRI_PER_SPHERE = raw(1, 1)
        raw = raw(2 : size(raw, 1));    % remove the triangle size part
        fclose(fid);
    
        for triangle_index = 1 : size(raw, 1)/9
            dp = raw((triangle_index - 1) * 9 + 1: triangle_index * 9);
            triangle.p1 = [dp(1:3)];
            triangle.p2 = [dp(4:6)];
            triangle.p3 = [dp(7:9)];
            
            if (triangle_index ==1 || mod(triangle_index - 1, TRI_PER_SPHERE) ~= 0)
                triangle_list = [triangle_list, triangle];
            else
                obj.triangle_list = triangle_list;
                obj_list = [obj_list, obj];
                triangle_list = [triangle];
            end
        end
        obj.triangle_list = triangle_list;
        obj_list = [obj_list, obj];   % wrap up
    else
        fprintf('File %s does not exist.\Nframe', fname);    % throw exception
    end
end

% for each frame, generate the sphere and render to the figure according to the extracted data
function perFrameRender(obj_list)
    objCount = size(obj_list, 2);

    for objId=1:objCount
        sphere = obj_list(objId).triangle_list;
        % extract sphere data point
        triCount = size(sphere, 2);

        % format the list of data points
        X = zeros(triCount, 3);
        Y = zeros(triCount, 3);
        Z = zeros(triCount, 3);
        
        % for each triangle inside the sphere
        for triId = 1:triCount
            triangle = sphere(triId);
            X(triId, 1) = triangle.p1(1);
            X(triId, 2) = triangle.p2(1);
            X(triId, 3) = triangle.p3(1);
            
            Y(triId, 1) = triangle.p1(2);
            Y(triId, 2) = triangle.p2(2);
            Y(triId, 3) = triangle.p3(2);
            
            Z(triId, 1) = triangle.p1(3);
            Z(triId, 2) = triangle.p2(3);
            Z(triId, 3) = triangle.p3(3);
        end
        
        % randomrize a color for the current sphere
        % C = rand(triCount, 3);

        % render triangle
        % fill3(X, Y, Z, C);
        color = 'r';
        if objId < objCount / 2
            color = 'b';
        end
        fill3(X, Y, Z, color);
        xlim([-0.6, 0.6]);
        ylim([-1, 1]);
        zlim([-1, 1]);
        hold on;
    end
end

% render or record the simulation video
function output(M, mode)
    if (mode == 1)
        % display the simulation movie
        movie(M, 3)
    else
        % record movie to avi
        v = VideoWriter('simulation.avi');
        v.FrameRate = 5;
        open(v);
        writeVideo(v,M);
        close(v);
    end
end