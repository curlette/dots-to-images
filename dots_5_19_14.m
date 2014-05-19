%Crops MOR1, GFP, and CD11 images in squares around tetrode tips. Uses
%cropped MOR1 images in WMB enhancement.
%
%INPUT: name(s) of folder(s) containing all relevant images for one
%hemisphere.
%Names of subfolders within hemisphere folders (one per folder) containing
%the dots, MOR1, GFP, etc. image folders
%Name(s) of corresponding Excel spreadsheets in which to save
%tetrode information for each hemisphere.
%save_format (format of output files)
%pixels_per_mm
%dim: side length of square used in cropping
%
%PRECONDITION: File structure [hemisphere folder] > [subfolder] > [folders
%for particular image types]
%
%OUTPUT: Folders of the format 'Output cropped images - tt[tetrode number]'
%containing subfolders MOR1, MOR1_0to99, MOR1_1to99, MOR1_0to0157, CD11_0to99, 
%Cd11_0to0157, CD11, GFP, WMB, and WM_BW

clear all
close all

foldernames = {'Dots input new'}; %fill in with names of 2 folders to process
subfoldernames = {'all aligned'}; %subfolders of hemisphere folders containing dots, MOR1, GFP, etc fodlers
tetrode_spreadsheet_names = {'tt_spreadsheet.xlsx'}; %fill in tetrode spreadsheet
%names for these 2 folders. Include '.xlsx' at end
pixels_per_mm = 985; %fill in from image
conversion_factor = pixels_per_mm/1000; %pixels per mm to pixels per micron
dim = 500*conversion_factor; %microns -> pixels
save_format = 'png';
ventricle_thresh = 3;
ventricle_diam = 100;

ventricle_diam = round(ventricle_diam*pixels_per_mm/1000);

for u = 1:length(foldernames)
    imagefolders = dir(fullfile([foldernames{u},'/',subfoldernames{u}]));
    
    dotsfolder = '';
    mor1folder = '';
    gfpfolder = '';
    cd11folder = '';
    
    %find image subfolders
    for k = 1:length(imagefolders)
        [path, name, ext] = fileparts(imagefolders(k).name);
        if isempty(strfind(imagefolders(k).name,'dots')) == 0
            dotsfolder = name;
        elseif isempty(strfind(imagefolders(k).name,'MOR1')) == 0
            mor1folder = name;
        elseif isempty(strfind(imagefolders(k).name,'GFP')) == 0
            gfpfolder = name;
        elseif isempty(strfind(imagefolders(k).name,'CD11')) == 0
            cd11folder = name;
        end
    end
    
    %cells to store different types of images and their names
    gfp = {};
    gfpnames = {};
    mor1 = {};
    mor1names = {};
    cd11 = {};
    cd11names = {};
    tt = {};
    ttnames = {};
    
    %finds and stores image files and names
    if strcmp(dotsfolder,'') == 0 %dots subfolder exists
        files = dir(fullfile([foldernames{u},'/',subfoldernames{u},'/',dotsfolder]));
        for k = 1:length(files)
            if isempty(strfind(files(k).name,'.tif')) == 0 || isempty(strfind(files(k).name,'.png')) == 0
                [path, name, ext] = fileparts(files(k).name);
                display(fullfile([foldernames{u},'/',subfoldernames{u},'/',dotsfolder,'/',files(k).name]))
                img = imread(fullfile([foldernames{u},'/',subfoldernames{u},'/',dotsfolder,'/',files(k).name]));
                tt{end + 1} = img;
                ttnames{1,end + 1} = name;
            end
        end
    end
    if strcmp(mor1folder,'') == 0
        files = dir(fullfile([foldernames{u},'/',subfoldernames{u},'/',mor1folder]));
        for k = 1:length(files)
            if isempty(strfind(files(k).name,'.tif')) == 0 || isempty(strfind(files(k).name,'.png')) == 0
                [path, name, ext] = fileparts(files(k).name);
                display(fullfile([foldernames{u},'/',subfoldernames{u},'/',mor1folder,'/',files(k).name]))
                img = imread(fullfile([foldernames{u},'/',subfoldernames{u},'/',mor1folder,'/',files(k).name]));
                mor1{end + 1} = img;
                mor1names{1,end + 1} = name;
            end
        end
    end
    if strcmp(gfpfolder,'') == 0
        files = dir(fullfile([foldernames{u},'/',subfoldernames{u},'/',gfpfolder]));
        for k = 1:length(files)
            if isempty(strfind(files(k).name,'.tif')) == 0 || isempty(strfind(files(k).name,'.png')) == 0
                [path, name, ext] = fileparts(files(k).name);
                display(fullfile([foldernames{u},'/',subfoldernames{u},'/',gfpfolder,'/',files(k).name]))
                img = imread(fullfile([foldernames{u},'/',subfoldernames{u},'/',gfpfolder,'/',files(k).name]));
                gfp{end + 1} = img;
                gfpnames{1,end + 1} = name;
            end
        end
    end
    if strcmp(cd11folder,'') == 0
        files = dir(fullfile([foldernames{u},'/',subfoldernames{u},'/',cd11folder]));
        for k = 1:length(files)
            if isempty(strfind(files(k).name,'.tif')) == 0 || isempty(strfind(files(k).name,'.png')) == 0
                [path, name, ext] = fileparts(files(k).name);
                display(fullfile([foldernames{u},'/',subfoldernames{u},'/',cd11folder,'/',files(k).name]))
                img = imread(fullfile([foldernames{u},'/',subfoldernames{u},'/',cd11folder,'/',files(k).name]));
                cd11{end + 1} = img;
                cd11names{1,end + 1} = name;
            end
        end
    end
    
    %matrices to store region numbers
    mor1nums = zeros(1,length(mor1));
    gfpnums = zeros(1,length(gfp));
    cd11nums = zeros(1,length(cd11));
    ttregnums = zeros(1,length(tt));
    ttregnums_no_OF = [];
    ttnums = cell(1,length(tt));
    ttnums_no_OF = {};
    
    %parse first 2 digits from image name string into region number for each
    %image type
    for k = 1:length(mor1nums)
        mor1nums(k) = str2num(mor1names{1,k}(1:2));
    end;
    for k = 1:length(gfpnums)
        gfpnums(k) = str2num(gfpnames{1,k}(1:2));
    end;
    for k = 1:length(cd11nums)
        cd11nums(k) = str2num(cd11names{1,k}(1:2));
    end;

    %store number from filename for tt images = tt number
    OFN = 1;
    for k = 1:length(ttnums)
        current = ttnames{1,k};
        index = strfind(current,'tt');
        if ~isempty(index)
            current_end = current(index:end);
            current_end_spaces = current_end == ' ';
            current_end_underscores = current_end == '_';
            current_end_either = current_end_spaces | current_end_underscores;
            current_end_either = find(current_end_either);
            if ~isempty(current_end_either)
                ttnum_end = current_end_either(1) - 1;
            else
                ttnum_end = length(current_end);
            end
            ttnums{k} = current_end(1:ttnum_end);
            ttnums_no_OF{end+1} = ttnums{k};
            ttregnums(k) = str2num(current(1:2));
            ttregnums_no_OF(end+1) = ttregnums(k);
        else
            ttnums{k} = ['OF',num2str(OFN)];
            ttregnums(k) = str2num(current(1:2));
            OFN = OFN + 1;
        end
    end;
    
    xlswrite(fullfile(foldernames{u},tetrode_spreadsheet_names{u}),ttnums_no_OF);
    xlswrite(fullfile(foldernames{u},tetrode_spreadsheet_names{u}),ttregnums_no_OF,1,'A2');
    
    %loop through all dot images
    for k = 1:length(ttregnums)
        mkdir(fullfile(foldernames{u},['Output cropped images - ',ttnums{k}]));
        mkdir(fullfile(foldernames{u},['Output cropped images - ',ttnums{k}],'MOR1'));
        mkdir(fullfile(foldernames{u},['Output cropped images - ',ttnums{k}],'GFP'));
        mkdir(fullfile(foldernames{u},['Output cropped images - ',ttnums{k}],'CD11'));
        
        mkdir(fullfile(foldernames{u},['Output cropped images - ',ttnums{k}],'MOR1_0to99'));
        mkdir(fullfile(foldernames{u},['Output cropped images - ',ttnums{k}],'CD11_0to99'));
        
        mkdir(fullfile(foldernames{u},['Output cropped images - ',ttnums{k}],'MOR1_1to99'));
        
        mkdir(fullfile(foldernames{u},['Output cropped images - ',ttnums{k}],'MOR1_0to0157'));
        mkdir(fullfile(foldernames{u},['Output cropped images - ',ttnums{k}],'CD11_0to0157'));
        
        mkdir(fullfile(foldernames{u},['Output cropped images - ',ttnums{k}],'WMB'));
        mkdir(fullfile(foldernames{u},['Output cropped images - ',ttnums{k}],'WM_BW'));
        
        img = tt{k};
        img = rgb2gray(img);
        img = img > 0; %~= 255; change based on whether dot image is light
        %on dark or dark on light
        [x,y,q] = size(img);
        measurements = regionprops(img,'centroid');
        [dot] = measurements.Centroid; %X,Y coordinates of dot centroid
        xpercent = dot(1)/y; %percentage of way along cols in image - so can
        %be mapped onto other images
        ypercent = dot(2)/x; %percentage of way along rows
        
        %loops through 3 layers above and below for cropping
        for j = ttregnums(k)-3:ttregnums(k)+3
            %as long as there's a mor1 image for region j, maps tetrode
            %coordinate to its dimensions and crops & saves [dim] x [dim] micron square
            if any(mor1nums == j)
                index = (mor1nums == j);
                index = find(index);
                mor1img = mor1{index};
                mor1name = mor1names{index};
                [a,b,q] = size(mor1img);
                tetrodex = round(ypercent*a);
                tetrodey = round(xpercent*b);
                newmor1 = mor1img((tetrodex - round(1/2*dim)):(tetrodex + ...
                    round(1/2*dim)),(tetrodey - round(1/2*dim)):(tetrodey + round(1/2*dim)),1:q);
                
                dark = newmor1 <= ventricle_thresh;
                ventricles = bwareaopen(dark, round((ventricle_diam/2)^2*pi));
                
                %image adjusted from 0% to 99%
                mor1dbl = im2double(newmor1);
                mor1novent = mor1dbl;
                mor1novent(ventricles) = [];
                mor1sort = sort(mor1novent);
                a = mor1sort(1);
                b = mor1sort(round(.99*length(mor1sort)));
                mor1_0to99 = imadjust(newmor1,[a;b],[0;1]);
                
                %image adjusted from 1% to 99%
                a = mor1sort(round(.01*length(mor1sort)));
                b = mor1sort(round(.99*length(mor1sort)));
                mor1_1to99 = imadjust(newmor1,[a;b],[0;1]);
                
                %BW image with 10px diameter white dot in center (tetrode
                %tip)
                if j == ttregnums(k)
                    bwdotimg = zeros(size(newmor1));
                    [a,b,q] = size(newmor1);
                    shapeInserter = vision.ShapeInserter('Shape','Circles','Fill',1,'FillColor','Custom','CustomFillColor',[255 255 255]);
                    circle = uint32([a/2, b/2, 5]);
                    %bwdotimg = repmat(double(bwdotimg)./255,[1 1 3]);
                    bwdotimg = step(shapeInserter, bwdotimg, circle);
                end
                
                %image adjusted from 0 to .157
                mor1_0to0157 =  imadjust(newmor1,[0;.157],[0;1]);
                
                imwrite(newmor1,fullfile(foldernames{u},['Output cropped images - ',ttnums{k}],'MOR1',[mor1name,'_',...
                    num2str(ttregnums(k)),'_',ttnums{k},'.',save_format]),save_format);
                imwrite(mor1_0to99,fullfile(foldernames{u},['Output cropped images - ',ttnums{k}],'MOR1_0to99',[mor1name,'_',...
                    num2str(ttregnums(k)),'_',ttnums{k},'.',save_format]),save_format);
                imwrite(mor1_1to99,fullfile(foldernames{u},['Output cropped images - ',ttnums{k}],'MOR1_1to99',[mor1name,'_',...
                    num2str(ttregnums(k)),'_',ttnums{k},'.',save_format]),save_format);
                imwrite(mor1_0to0157,fullfile(foldernames{u},['Output cropped images - ',ttnums{k}],'MOR1_0to0157',[mor1name,'_',...
                    num2str(ttregnums(k)),'_',ttnums{k},'.',save_format]),save_format);
                if j == ttregnums(k)
                    imwrite(bwdotimg,fullfile(foldernames{u},['Output cropped images - ',ttnums{k}],[mor1name,'_',...
                        num2str(ttregnums(k)),'_',ttnums{k},'_bwdotimg','.',save_format]),save_format);
                end
                
            end
            
            if any(gfpnums == j)
                index = (gfpnums == j);
                index = find(index);
                gfpimg = gfp{index};
                gfpname = gfpnames{index};
                [c,d,s] = size(gfpimg);
                tetrodex = round(ypercent*c);
                tetrodey = round(xpercent*d);
                newgfp = gfpimg((tetrodex - round(1/2*dim)):(tetrodex + ...
                    round(1/2*dim)),(tetrodey - round(1/2*dim)):(tetrodey + round(1/2*dim)),1:q);
                imwrite(newgfp,fullfile(foldernames{u},['Output cropped images - ',ttnums{k}],'GFP',[gfpname,'_',...
                    num2str(ttregnums(k)),'_',ttnums{k},'.',save_format]),save_format);
            end;
            if any(cd11nums == j)
                index = (cd11nums == j);
                index = find(index);
                cd11img = cd11{index};
                cd11name = cd11names{index};
                [e,f,r] = size(cd11img);
                tetrodex = round(ypercent*e);
                tetrodey = round(xpercent*f);
                size(cd11img)
                newcd11 = cd11img((tetrodex - round(1/2*dim)):(tetrodex + ...
                    round(1/2*dim)),(tetrodey - round(1/2*dim)):(tetrodey + round(1/2*dim)),1:q);
                
                cd11dbl = im2double(newcd11);
                
                cd11novent = cd11dbl;
                cd11novent(ventricles) = [];
                cd11sort = sort(cd11novent);
                a = cd11sort(1);
                b = cd11sort(round(.99*length(cd11sort)));
                cd11_0to99 = imadjust(newcd11,[a;b],[0;1]);
                
                cd11_0to0157 =  imadjust(newcd11,[0;.157],[0;1]);
                
                imwrite(newcd11,fullfile(foldernames{u},['Output cropped images - ',ttnums{k}],'CD11',[cd11name,'_',...
                    num2str(ttregnums(k)),'_',ttnums{k},'.',save_format]),save_format);
                imwrite(cd11_0to99,fullfile(foldernames{u},['Output cropped images - ',ttnums{k}],'CD11_0to99',[cd11name,'_',...
                    num2str(ttregnums(k)),'_',ttnums{k},'.',save_format]),save_format);
                imwrite(cd11_0to0157,fullfile(foldernames{u},['Output cropped images - ',ttnums{k}],'CD11_0to0157',[cd11name,'_',...
                    num2str(ttregnums(k)),'_',ttnums{k},'.',save_format]),save_format);
            end;
        end;
    end;
    
    %---- white matter enhancement ----
    low_in_increment = 0.001;
    high_in_increment = 0.001;
    nonblackfrac_target = 0.15;
    nonblack_areathresh_percent = 0.1;
    medPixVal_target = 0.5;
    
    for k = 1:length(ttregnums)
        files1 = dir(fullfile(foldernames{u},['Output cropped images - ',ttnums{k},'/MOR1'],['*.',save_format]));
        num_files = length(files1);
        
        images = {}; %cell to store images
        image_names = {}; %cell to store image names
        
        for r = 1:num_files
            filename = files1(r).name;
            [pathstr, name, ext] = fileparts(filename);
            img = imread(fullfile(foldernames{u},['Output cropped images - ',ttnums{k},'/MOR1'],filename));
            images{end + 1} = im2double(img);
            image_names{end + 1} = name;
        end;
        
        for m = 1:length(images)
            current = images{m};
            %current = rgb2gray(current);
            [a,b,q] = size(current);
            
            black = current == 0;
            
            pixlists = regionprops(black, 'PixelList');
            [numRegions, ~] = size(pixlists);
            
            blackedge = {};
            
            %finds black borders
            for p = 1:numRegions
                edge = 0;
                currentPixList = pixlists(p).PixelList;
                [currentUnwantedArea, ~] = size(currentPixList);
                for n = 1:currentUnwantedArea
                    curY = currentPixList(n,1);
                    curX = currentPixList(n,2);
                    if curY == 1 || curY == b || curX == 1 || curX == a
                        edge = 1;
                    end
                end
                if edge
                    blackedge{end + 1} = currentPixList;
                end
            end
            
            currentinv = imcomplement(current);
            
            numignoredpx_edges = 0;
            edge_array = zeros(a,b);
            for p = 1:length(blackedge)
                currentPixList = blackedge{p};
                [currentUnwantedArea, ~] = size(currentPixList);
                numignoredpx_edges = numignoredpx_edges + currentUnwantedArea;
                for n = 1:currentUnwantedArea
                    curY = currentPixList(n,1);
                    curX = currentPixList(n,2);
                    currentinv(curX,curY) = 0;
                    edge_array(curX,curY) = 1;
                end
            end
            newname5 = [image_names{m},'_edge.',save_format];
            imwrite(edge_array,fullfile(foldernames{u},['Output cropped images - ',ttnums{k}],newname5),save_format);
            
            low_in = 1;
            high_in = 1;
            nonblackfrac = 0;
            nonblack_areathresh = nonblack_areathresh_percent*a*b;
            current2 = zeros(a,b);
            nonblackimg = zeros(a,b);
            largewm = zeros(a,b);
            nonblackimg2 = zeros(a,b);
            
            %decreases low_in until nonblackfrac >= nonblackfrac_target and makes large WM/ventricles black in images to be saved
            while nonblackfrac < nonblackfrac_target
                low_in_previous = low_in;
                low_in = low_in - low_in_increment;
                if low_in < 0
                    low_in = 0;
                    break
                end
                current2_previous = current2;
                current2 = imadjust(currentinv,[low_in;high_in], [0;1]);
                nonblackimg_previous = nonblackimg;
                nonblackimg = current2 > 0;
                pixlists = regionprops(nonblackimg, 'PixelList');
                [numRegions, ~] = size(pixlists);
                numignoredpx2 = numignoredpx_edges;
                largewm_previous = largewm;
                largewm = zeros(a,b);
                for p = 1:numRegions
                    edge = 0;
                    currentPixList = pixlists(p).PixelList;
                    [currentUnwantedArea, ~] = size(currentPixList);
                    if currentUnwantedArea > nonblack_areathresh
                        for n = 1:currentUnwantedArea
                            curY = currentPixList(n,1);
                            curX = currentPixList(n,2);
                            if curY == 1 || curY == b || curX == 1 || curX == a
                                edge = 1;
                            end
                        end
                        if edge
                            numignoredpx2 = numignoredpx2 + currentUnwantedArea;
                            for n = 1:currentUnwantedArea
                                curY = currentPixList(n,1);
                                curX = currentPixList(n,2);
                                largewm(curX,curY) = 1;
                            end
                        end
                    end
                end
                nonblackimg = current2 > 0;
                nonblackimg2_previous = nonblackimg2;
                nonblackimg2 = nonblackimg;
                nonblackimg2(logical(largewm)) = 0;
                nonblackfrac_previous = nonblackfrac;
                nonblackfrac = sum(nonblackimg2(:))/(a*b - numignoredpx2);
            end
            if abs(nonblackfrac_target - nonblackfrac_previous) < abs(nonblackfrac_target - nonblackfrac)
                current2 = current2_previous;
                low_in = low_in_previous;
                nonblackimg = nonblackimg_previous;
                largewm = largewm_previous;
                nonblackfrac = nonblackfrac_previous;
                nonblackimg2 = nonblackimg2_previous;
            end
            
            %saves all wm image before large regions are removed
            newname3 = [image_names{m},'_wm_bw.',save_format];
            imwrite(nonblackimg2,[foldernames{u},'/Output cropped images - ',ttnums{k},'/WM_BW/',newname3], save_format);
            
            %decreases high_in until medPixVal >= medPixVal_target
            medPixVal = median(current2(logical(nonblackimg2)));
            if medPixVal < medPixVal_target
                while medPixVal < medPixVal_target
                    high_in_previous = high_in;
                    high_in = high_in - high_in_increment;
                    if high_in <= low_in
                        high_in = high_in_previous;
                        break
                    end
                    current2_previous = current2;
                    current2 = imadjust(currentinv,[low_in;high_in], [0;1]);
                    current2(logical(largewm)) = 0;
                    medPixVal_previous = medPixVal;
                    medPixVal = median(current2(nonblackimg2));
                end
                if abs(medPixVal_target - medPixVal_previous) < abs(medPixVal_target - medPixVal)
                    medPixVal = medPixVal_previous;
                    high_in = high_in_previous;
                    current2 = current2_previous;
                end
            end
            figure; imshow(current2);
            %saves output images
            newname4 = [image_names{m},'_large_wm_bw.',save_format];
            imwrite(largewm,[foldernames{u},'/Output cropped images - ',ttnums{k},'/WM_BW/',newname4], save_format);
            
            newname = [image_names{m},'_wmb.',save_format];
            newname2 = [image_names{m},'_small_wm_bw.',save_format];
            imwrite(current2,[foldernames{u},'/Output cropped images - ',ttnums{k},'/WMB/',newname], save_format);
            imwrite(nonblackimg,[foldernames{u},'/Output cropped images - ',ttnums{k},'/WM_BW/',newname2], save_format);
            
            display(nonblackfrac);
            display(medPixVal);
        end
    end
end