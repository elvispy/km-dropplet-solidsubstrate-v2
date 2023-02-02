function export_animation(varargin)
    if nargin == 1
       
    else
        [file, path] = uigetfile("*.mat");
        fullfilepath = fullfile(path, file);
        load(fullfilepath, "recorded_conditions", "length_unit", "velocity_unit", "presure_unit");
        vidObj = VideoWriter(replace(fullfilepath, ".mat", ".mp4"), "MPEG-4");
        set(vidObj, 'Quality', 100, 'FrameRate', 100);
        open(vidObj);
        
        for ii = 1:size(recorded_conditions, 1)
            recorded_conditions{ii}.amplitude_defor = 1;
            plot_condition(1, recorded_conditions{ii}, 5);
            writeVideo(vidObj, getframe(gcf));
            

        end
        close(vidObj);
    end

end