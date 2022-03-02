function [D_3d] = fx_Scrollable_3d_view(varargin)
% Graphical User Interface for Scrolling through a 3D Image, in 3 dimension 
% slice.
%% 
% [D_3d] = Scrollable_3d_view(image,Dim)
% gets a 3d image as input and you can initialize which slice you want to 
% start with by giving a 1 in 3 vector to the function as "Dim" input 
% input argument.

% [D_3d] = Scrollable_3d_view(image)
% gets a 3d image as input and initialize with the slice equal to half of
% the dimention.

% As an output (D_3d). it gives you the slice number for the last slices
% that you were looking at. You can use the "Done" buttom to get the output.

% Author: Hamed Hooshangnejad
% E-mail: hh1368hh@gmail.com
% Release: 1.0
% Release date: 6/18/19



if size(varargin,2)==1
    File=varargin{1};
    D_3d=round(size(File)/2);
else
    File=varargin{1};
    D_3d=varargin{2};
end


global fdata loaded ax img img_rs img_rs_slice d
img_rs=File;
d=D_3d;

fdata = figure( 'windowstyle', 'normal', 'resize', 'on', 'position',[200 100 900 600],...
    'name', 'Visualizer', 'NumberTitle', 'off','KeyPressFcn',@keyPress);

loadfile;

done = uicontrol('parent', fdata, 'Style', 'pushbutton', 'units', 'normalized', 'String', 'Done', ...
    'Position', [0.03 0.8 0.05 0.05], 'Callback', @done_loading);


    function done_loading(hObject, eventdata)
        D_3d=d;

        close(fdata)
        
        uiresume;
    end


uiwait;
end

function loadfile(hObject, eventdata)
global fdata img loaded ax d info Slider_a Slider_s Slider_c text_a text_s text_c img_rs img_rs_slice s_h I_S I_C I_A
global dline_cell vis_conts colors
%% Loading the Image

img_axial=img_rs(:,:,d(3));
img_sagitall=(squeeze(img_rs(:,d(2),:)))';
img_coronal=(squeeze(img_rs(d(1),:,:)))';
d_line_3d=cell2mat(dline_cell);


s_h=slice(double(img_rs),d(1),d(2),d(3));
% axis image

% changed by yxy
mymap = jet;
% mymap(1, :) = [0 0 0]; % set the minmun to black
colormap(mymap);
colorbar;

shading interp;
set(gca,'Zdir','reverse')
view(-150,35)
if ~isempty(d_line_3d)
    hold on
    h3d=plot3(d_line_3d(:,1),d_line_3d(:,2),d_line_3d(:,3), 'LineStyle', 'none','Marker','.','Color','r');
    hold off
end
axis on;
% added by yxy
xlabel('x');
ylabel('y');
zlabel('z');

if ~isempty(vis_conts)
    for con_k=1:length(vis_conts)
        BW=vis_conts{con_k}.mask_rs;
        figure(fdata)
        coor=find(BW, 1);
        
        if ~isempty(coor)

            con_B(con_k)=patch(isosurface(BW,0),'FaceColor',colors(con_k,:),'EdgeColor','none');
            lighting phong;
            camlight;
            hold off
            
        end
    end
    
end
%% creating the sliders

Slider_a = uicontrol('parent', fdata, 'Style', 'slider', 'units', 'normalized',...
    'Position', [0.13 0.05 0.2 0.03],'Callback', @slider_a_func);
set(Slider_a,'Sliderstep',[1/(size(img_rs,3)-1) 10/(size(img_rs,3)-1)],...
            'max',size(img_rs,3),'min',1,'Value',d(3))
text_a= uicontrol('parent', fdata, 'Style', 'text', 'String', num2str(get(Slider_a,'value')), 'units', 'normalized', ...
    'FontSize',7,'Position', [0.2 0.0 0.05 0.05]);

text_a1= uicontrol('parent', fdata, 'Style', 'text', 'String','Z-Axis', 'units', 'normalized', ...
    'FontSize',7,'Position', [0.2 0.08 0.05 0.02]);
%==================================
Slider_s = uicontrol('parent', fdata, 'Style', 'slider', 'units', 'normalized',...
    'Position', [0.4 0.05 0.2 0.03],'Callback', @slider_s_func);
set(Slider_s,'Sliderstep',[1/(size(img_rs,2)-1) 10/(size(img_rs,2)-1)],...
            'max',size(img_rs,2),'min',1,'Value',d(2))
text_s= uicontrol('parent', fdata, 'Style', 'text', 'String', num2str(get(Slider_s,'value')), 'units', 'normalized', ...
    'FontSize',7,'Position', [0.48 0.0 0.05 0.05]);
text_s1= uicontrol('parent', fdata, 'Style', 'text', 'String','X-Axis', 'units', 'normalized', ...
    'FontSize',7,'Position', [0.48 0.08 0.05 0.02]);
%==================================
Slider_c = uicontrol('parent', fdata, 'Style', 'slider', 'units', 'normalized',...
    'Position', [0.7 0.05 0.2 0.03],'Callback', @slider_c_func);
set(Slider_c,'Sliderstep',[1/(size(img_rs,1)-1) 10/(size(img_rs,1)-1)],...
            'max',size(img_rs,1),'min',1,'Value',d(1))
text_c= uicontrol('parent', fdata, 'Style', 'text', 'String', num2str(get(Slider_c,'value')), 'units', 'normalized', ...
    'FontSize',7,'Position', [0.78 0.0 0.05 0.05]);
text_c1= uicontrol('parent', fdata, 'Style', 'text', 'String','Y-Axis', 'units', 'normalized', ...
    'FontSize',7,'Position', [0.78 0.08 0.05 0.02]);
loaded=1;
end

function slider_a_func(hObject, eventdata)
global fdata img loaded ax d info Slider_a Slider_s Slider_c text_a text_s text_c img_rs s_h   I_A

% new value of d for visualizing
nd=get(Slider_a,'value');
set(text_a, 'String', num2str(round(nd)));

img_axial=img_rs(:,:,round(nd));

d(3)=round(nd);

s_h(3).CData=img_axial;
newZ=round(nd)*ones(size(s_h(3).ZData));
s_h(3).ZData=newZ;

end

function slider_s_func(hObject, eventdata)
global fdata img loaded ax d info Slider_a Slider_s Slider_c text_a text_s text_c  pixel_spc slide_spc img_rs s_h I_S

% new value of d for visualizing
nd=get(Slider_s,'value');
set(text_s, 'String', num2str(round(nd)));
img_sagitall=(squeeze(img_rs(:,round(nd),:)))';


d(2)=round(nd);


s_h(1).CData=img_sagitall';
newX=round(nd)*ones(size(s_h(1).XData));
s_h(1).XData=newX;

end

function slider_c_func(hObject, eventdata)
global fdata img loaded ax d info Slider_a Slider_s Slider_c text_a text_s text_c  pixel_spc slide_spc img_rs s_h I_C

% new value of d for visualizing
nd=get(Slider_c,'value');
set(text_c, 'String', num2str(round(nd)));
img_coronal=(squeeze(img_rs(round(nd),:,:)))';

d(1)=round(nd);

s_h(2).CData=img_coronal';
newY=round(nd)*ones(size(s_h(2).YData));
s_h(2).YData=newY;
end
