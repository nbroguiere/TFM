function stackshow(stack)
% Displays 3D stacks with a slider, can be binary or 2 or three colors. 
% Adapted from https://stackoverflow.com/questions/28256106/image-stack-display-in-matlab-using-a-slider

stack=squeeze(stack);
s=size(stack);
if length(s)==4
    if s(3)==3
        disp(['Color stack detected with ' num2str(s(4)) ' layers (color is third dim, layer 4th dim).']);
    elseif s(3)==2
        disp('Only two color channels, interpret as red and green');
        stack=double(stack);
        stack=stack./max(stack(:));
        stack=cat(3,stack,zeros(s(1),s(2),1,s(4)));
    elseif s(4)==3
        stack=permute(stack,[1,2,4,3]);
        stack=double(stack);
        stack=stack./max(stack(:));
        disp(['Color stack detected with ' num2str(s(3)) ' layers (color is 4th dim, layer 3rd dim).']);
    elseif s(4)==2
        stack=permute(stack,[1,2,4,3]);
        stack=double(stack);
        stack=stack./max(stack(:));
        stack=cat(3,stack,zeros(s(1),s(2),1,s(3)));
        disp(['Only two color channels, interpret as red and green. ' num2str(s(3)) ' layers (color is 4th dim, layer 3rd dim).']);
    else
        disp('Dimensions not supported.');
        return
    end
elseif length(s)==3
    if islogical(stack)
        disp('Black and white stack detected')
        stack=reshape(stack,[s(1) s(2) 1 s(3)]);
        stack=double(cat(3,stack,stack,stack));
    elseif strcmp(class(stack),'uint8') || strcmp(class(stack),'uint16') 
        disp('Integer gray level image detected, normalizing to highest value')
        stack=double(stack);
        stack=stack/max(stack(:));
        stack=reshape(stack,[s(1) s(2) 1 s(3)]);
        stack=cat(3,stack,stack,stack);
    elseif strcmp(class(stack),'double')
        disp('Double gray level image detected, normalizing to highest value')
        stack=(stack-min(stack(:)))/(max(stack(:))-min(stack(:)));
        stack=reshape(stack,[s(1) s(2) 1 s(3)]);
        stack=cat(3,stack,stack,stack);
    end
elseif length(s)==2
    disp('detected single plane picture');
    imshow(stack);
    return
else
    disp('Stack has uncompatible dimensions');
    return
end
NumFrames = size(stack,4); %// Check below for dummy 4D matrix/image sequence
screen_size = [100 100 1700 900]; 
%hFig = figure('Position',[100 100 500 500],'Units','normalized');
hFig = figure('Position',screen_size,'Units','normalized');
handles.axes1 = axes('Units','normalized','Position',[.2 .2 .6 .6]);

%// Create slider and listener object for smooth visualization
handles.SliderFrame = uicontrol('Style','slider','Position',[0 0 1700 20],'Min',1,'Max',NumFrames,'Value',1,'SliderStep',[1/NumFrames 2/NumFrames],'Callback',@XSliderCallback);
handles.SliderxListener = addlistener(handles.SliderFrame,'Value','PostSet',@(s,e) XListenerCallBack);

%handles.Text1 = uicontrol('Style','Text','Position',[180 420 60 30],'String','Current frame');
handles.Edit1 = uicontrol('Style','Edit','Position',[0 20 30 30],'String','1');

%// Use setappdata to store the image stack and in callbacks, use getappdata to retrieve it and use it. Check the docs for the calling syntax.

setappdata(hFig,'MyMatrix',stack); %// You could use %//setappdata(0,'MyMatrix',MyMatrix) to store in the base workspace. 

%// Display 1st frame
imshow(stack(:,:,:,1))

%// IMPORTANT. Update handles structure.
guidata(hFig,handles);

%// Listener callback, executed when you drag the slider.

    function XListenerCallBack

        %// Retrieve handles structure. Used to let MATLAB recognize the
        %// edit box, slider and all UI components.
        handles = guidata(gcf);

%// Here retrieve MyMatrix using getappdata.
stack = getappdata(hFig,'MyMatrix');

        %// Get current frame
        CurrentFrame = round((get(handles.SliderFrame,'Value')));
        set(handles.Edit1,'String',num2str(CurrentFrame));

        %// Display appropriate frame.
        imshow(stack(:,:,:,CurrentFrame),'Parent',handles.axes1);

        guidata(hFig,handles);
    end


%// Slider callback; executed when the slider is release or you press
%// the arrows.
    function XSliderCallback(~,~)

        handles = guidata(gcf);

%// Here retrieve MyMatrix using getappdata.
    stack = getappdata(hFig,'MyMatrix');

        CurrentFrame = round((get(handles.SliderFrame,'Value')));
        set(handles.Edit1,'String',num2str(CurrentFrame));

        imshow(stack(:,:,:,CurrentFrame),'Parent',handles.axes1);

        guidata(hFig,handles);
    end

end