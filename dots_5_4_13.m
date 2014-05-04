%Crops MOR1, GFP, and CD11 images in squares around tetrode tips. Uses 
%cropped MOR1 images in WMB enhancement. 
%
%INPUT: name(s) of folder(s) containing all relevant images for one
%hemisphere. 
%Name(s) of corresponding Excel spreadsheets in which to save
%tetrode information for each hemisphere.
%save_format (format of output files)
%pixels_per_mm
%dim: side length of square used in cropping 
%
%PRECONDITION: All image types (MOR1,CD11, GFP, and dot files for each
%hemisphere) are in the same folder 
%
%OUTPUT: Folders of the format 'Output cropped images - tt[tetrode number]'
%containing subfolders MOR1, MOR1_0to99, MOR1_0to157, CD11_0to99, Cd11_0to157,
%CD11, GFP, WMB, and WM_BW

foldernames = {'all aligned'}; %fill in with names of 2 folders to process
tetrode_spreadsheet_names = {'tt_spreadsheet.xlsx'}; %fill in tetrode spreadsheet 
%names for these 2 folders. Include '.xlsx' at end

for u = 1:length(foldernames)
    tif_files = dir(fullfile(foldernames{u},'*.tif')); %tif files in folder
    png_files = dir(fullfile(foldernames{u},'*.png'));

    %cells to store different types of images and their names 
    gfp = {};
    gfpnames = {};
    mor1 = {};
    mor1names = {};
    cd11 = {};
    cd11names = {};
    tt = {};
    ttnames = {};

    pixels_per_mm = 1000; %fill in from image
    conversion_factor = pixels_per_mm/1000; %pixels per mm to pixels per micron
    dim = 1000*conversion_factor; %microns -> pixels
    save_format = 'png'; 

    %go through images and store mor1, gfp, cd11, and tt and their names
    for k = 1:length(tif_files)
        filename = tif_files(k).name;
        [pathstr, name, ext] = fileparts(filename);
        
        display(fullfile(foldernames{u},filename))
        img = imread(fullfile(foldernames{u},filename));

        if isempty(strfind(name, 'MOR1')) == 0 || isempty(strfind(name, '_0001_')) == 0
            mor1{end + 1} = img;
            mor1names{1,end + 1} = name;
        elseif isempty(strfind(name, 'GFP')) == 0 || isempty(strfind(name, '_0002_')) == 0
            gfp{end + 1} = img;
            gfpnames{1,end + 1} = name;
        elseif isempty(strfind(name, 'CD11')) == 0 || isempty(strfind(name, '_0000_')) == 0
            cd11{end + 1} = img;
            cd11names{1,end + 1} = name;
        elseif isempty(strfind(name, 'tt')) == 0
            tt{end + 1} = img;
            ttnames{1,end + 1} = name;
        end;
    end;

    for k = 1:length(png_files)
        filename = png_files(k).name;
        [pathstr, name, ext] = fileparts(filename);
        
        display(fullfile(foldernames{u},filename))
        img = imread(fullfile(foldernames{u},filename));

        if isempty(strfind(name, 'MOR1')) == 0 || isempty(strfind(name, '_0001_')) == 0
            mor1{end + 1} = img;
            mor1names{1,end + 1} = name;
        elseif isempty(strfind(name, 'GFP')) == 0 || isempty(strfind(name, '_0002_')) == 0
            gfp{end + 1} = img;
            gfpnames{1,end + 1} = name;
        elseif isempty(strfind(name, 'CD11')) == 0 || isempty(strfind(name, '_0000_')) == 0
            cd11{end + 1} = img;
            cd11names{1,end + 1} = name;
        elseif isempty(strfind(name, 'tt')) == 0
            tt{end + 1} = img;
            ttnames{1,end + 1} = name;
        end;
    end;

    %matrices to store region numbers
    mor1nums = zeros(1,length(mor1));
    gfpnums = zeros(1,length(gfp));
    cd11nums = zeros(1,length(cd11));
    ttregnums = zeros(1,length(tt));
    ttnums = cell(1,length(tt));

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
    for k = 1:length(ttregnums)
        ttregnums(k) = str2num(ttnames{1,k}(1:2));
    end;
    %store last number from filename for tt images = tt number
    for k = 1:length(ttnums)
        current = ttnames{1,k};
        index = strfind(current,'tt');
        ttnums{k} = current(index + 2:length(current)); %whatever's after the second t in 'tt'
    end;

    xlswrite(tetrode_spreadsheet_names{u},ttnums);
    xlswrite(tetrode_spreadsheet_names{u},ttregnums,1,'A2');

    %loop through all dot images 
    for k = 1:length(ttregnums)
        mkdir(['Output cropped images - tt',ttnums{k}]);
        mkdir(['Output cropped images - tt',ttnums{k},'/MOR1']);
        mkdir(['Output cropped images - tt',ttnums{k},'/GFP']);
        mkdir(['Output cropped images - tt',ttnums{k},'/CD11']);
        
        mkdir(['Output cropped images - tt',ttnums{k},'/MOR1_0to99']);       
        mkdir(['Output cropped images - tt',ttnums{k},'/CD11_0to99']);
        
        mkdir(['Output cropped images - tt',ttnums{k},'/MOR1_0to157']);        
        mkdir(['Output cropped images - tt',ttnums{k},'/CD11_0to157']);
        
        mkdir(['Output cropped images - tt',ttnums{k},'/WMB']);
        mkdir(['Output cropped images - tt',ttnums{k},'/WM_BW']);
        
        img = tt{k};
%         img = rgb2gray(img);
        img = img > 0; %~= 255; change based on whether dot image is light 
        %on dark or dark on light
        [x,y,q] = size(img); 
        measurements = regionprops(img,'centroid'); 
        [dot] = measurements.Centroid; %X,Y coordinates of dot centroid
        xpercent = dot(1)/y; %percentage of way along cols in image - so can
        %be mapped onto other images
        ypercent = dot(2)/x; %percentage of way along rows 

        %loops through 3 layers above and below for croppping
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
                
                %image adjusted from 0% to 99%
                mor1dbl = im2double(newmor1);
                mor1sort = sort(mor1dbl(:));
                a = mor1sort(1);
                b = mor1sort(round(.99*length(mor1sort)));
                mor1_0to99 = imadjust(newmor1,[a;b],[0;1]);   
                
                %image adjusted from 0 to .157
                mor1_0to157 =  imadjust(newmor1,[0;.157],[0;1]);                
                
                imwrite(newmor1,['Output cropped images - tt',ttnums{k},'/MOR1/',mor1name,'_',...
                    num2str(ttregnums(k)),'_tt',ttnums{k},'.tif'],save_format);
                imwrite(mor1_0to99,['Output cropped images - tt',ttnums{k},'/MOR1_0to99/',mor1name,'_',...
                    num2str(ttregnums(k)),'_tt',ttnums{k},'.tif'],save_format);
                imwrite(bwdotimg,['Output cropped images - tt',ttnums{k},'/MOR1_0to99/',mor1name,'_',...
                    num2str(ttregnums(k)),'_tt',ttnums{k},'_bwdotimg','.tif'],save_format);
                imwrite(mor1_0to157,['Output cropped images - tt',ttnums{k},'/MOR1_0to157/',mor1name,'_',...
                    num2str(ttregnums(k)),'_tt',ttnums{k},'.tif'],save_format);
            end;
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
                imwrite(newgfp,['Output cropped images - tt',ttnums{k},'/GFP/',gfpname,'_',...
                    num2str(ttregnums(k)),'_tt',ttnums{k},'.tif'],save_format);
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
                cd11sort = sort(cd11dbl(:));
                a = cd11sort(1);
                b = cd11sort(round(.99*length(cd11sort)));
                cd11_0to99 = imadjust(newcd11,[a;b],[0;1]);   
                
                cd11_0to157 =  imadjust(newcd11,[0;.157],[0;1]);
                
                imwrite(newcd11,['Output cropped images - tt',ttnums{k},'/CD11/',cd11name,'_',...
                    num2str(ttregnums(k)),'_tt',ttnums{k},'.tif'],save_format);
                imwrite(cd11_0to99,['Output cropped images - tt',ttnums{k},'/CD11_0to99/',cd11name,'_',...
                    num2str(ttregnums(k)),'_tt',ttnums{k},'.tif'],save_format);
                imwrite(cd11_0to157,['Output cropped images - tt',ttnums{k},'/CD11_0to157/',cd11name,'_',...
                    num2str(ttregnums(k)),'_tt',ttnums{k},'.tif'],save_format);
            end;
        end;
        %BW image with 10px diameter white dot in center (tetrode
        %tip)
        bwdotimg = zeros(size(newmor1));
        [a,b,q] = size(newmor1);
        shapeInserter = vision.ShapeInserter('Shape','Circles','Fill',1,'FillColor','Custom','CustomFillColor',[255 255 255]);            
        circle = uint32([a/2, b/2, 5]);
        bwdotimg = repmat(double(bwdotimg)./255,[1 1 3]);            
        bwdotimg = step(shapeInserter, bwdotimg, circle);
    end;

    %---- white matter enhancement ----
    N = 20;
    low_in_increment = 0.001;
    high_in_increment = 0.001;
    nonblackfrac_target = 0.2;
    nonblack_areathresh_percent = 0.1;
    medPixVal_target = 0.5;
    
    for k = 1:length(ttregnums)
        files1 = dir(fullfile(['Output cropped images - tt',ttnums{k},'/MOR1'],'*.tif')); 
        num_files = length(files1); 
        
        images = {}; %cell to store images
        image_names = {}; %cell to store image names

        for r = 1 : num_files
            filename = files1(r).name;
            [pathstr, name, ext] = fileparts(filename);            
            img = imread(fullfile(['Output cropped images - tt',ttnums{k},'/MOR1'],filename));            
            images{end + 1} = im2double(img);
            image_names{end + 1} = name;         
        end;

        for m = 1:length(images)
            current = images{m};

            [a,b] = size(current);
    
            black = current == 0; 

            pixlists = regionprops(black, 'PixelList');
            [numRegions, ~] = size(pixlists);

            blackedge = {};

            %finding black borders
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
                if edge == 1
                    blackedge{end + 1} = currentPixList;
                end
            end

            currentinv = imcomplement(current);

            numignoredpx = 0;

            for p = 1:length(blackedge)
                currentPixList = blackedge{p};
                [currentUnwantedArea, ~] = size(currentPixList);
                numignoredpx = numignoredpx + currentUnwantedArea;
                for n = 1:currentUnwantedArea
                    curY = currentPixList(n,1);
                    curX = currentPixList(n,2);
                    currentinv(curX,curY) = 0;
                end
            end

            low_in = 1;
            high_in = 1;
            nonblackfrac = 0;
            nonblack_areathresh = 10000;%nonblack_areathresh_percent*(a*b - numignoredpx);

            %decreases low_in until nonblackfrac >= nonblackfrac_target
            while nonblackfrac < nonblackfrac_target
                low_in = low_in - low_in_increment;
                current2 = imadjust(currentinv,[low_in;high_in], [0;1]);
                nonblackimg = current2 > 0;
                pixlists = regionprops(nonblackimg, 'PixelList');
                [numRegions, ~] = size(pixlists);
                numignoredpx2 = 0;
                for p = 1:numRegions
                    currentPixList = pixlists(p).PixelList;
                    [currentUnwantedArea, ~] = size(currentPixList);

                    if currentUnwantedArea > nonblack_areathresh
                       numignoredpx = numignoredpx + currentUnwantedArea;                
                    end            
                end
                nonblackimg = current2 > 0;
                low_in_previous = low_in;
                nonblackfrac_previous = nonblackfrac;
                nonblackfrac = sum(nonblackimg(:))/(a*b - numignoredpx2);
            end
            if abs(nonblackfrac_target - nonblackfrac_previous) < abs(nonblackfrac_target - nonblackfrac)
                nonblackfrac = nonblackfrac_previous;
                low_in = low_in_previous;
            end

            newname3 = [image_names{m},'_wm_bw.png'];
            imwrite(nonblackimg,['Output cropped images - tt',ttnums{k},'/WM_BW/',newname3], 'png');

            %make large WM/ventricles black in both images to be saved
            largewm = zeros(size(current2));
            nonblackimg = current2 > 0;
            pixlists = regionprops(nonblackimg, 'PixelList');
            [numRegions, ~] = size(pixlists);
            for p = 1:numRegions
                currentPixList = pixlists(p).PixelList;
                [currentUnwantedArea, ~] = size(currentPixList);

                if currentUnwantedArea > nonblack_areathresh            
                    for n = 1:currentUnwantedArea
                        curY = currentPixList(n,1);
                        curX = currentPixList(n,2);
                        current2(curX, curY) = 0;
                        nonblackimg(curX, curY) = 0;   
                        largewm(curX,curY) = 1;
                    end       
                end        
            end

            medPixVal = median(current2(current2 > 0));
            if medPixVal < medPixVal_target
                while medPixVal < medPixVal_target
                    high_in = high_in - high_in_increment;
                    if low_in >= high_in
                        break
                    end
                    current2 = imadjust(currentinv,[low_in;high_in], [0;1]);
                    high_in_previous = high_in;
                    medPixVal_previous = medPixVal;
                    medPixVal = median(current2(current2 > 0));
                end
                if abs(medPixVal_target - medPixVal_previous) < abs(medPixVal_target - medPixVal)
                    medPixVal = medPixVal_previous;
                    high_in = high_in_previous;
                end
            end   

            newname4 = [image_names{m},'_large_wm_bw.png'];
            imwrite(largewm,['Output cropped images - tt',ttnums{k},'/WM_BW/',newname4], 'png');

            newname = [image_names{m},'_wmb.png'];
            newname2 = [image_names{m},'_small_wm_bw.png'];   
            imwrite(current2,['Output cropped images - tt',ttnums{k},'/WMB/',newname], 'png');  
            imwrite(nonblackimg,['Output cropped images - tt',ttnums{k},'/WM_BW/',newname2], 'png');  

            display(nonblackfrac);
            display(medPixVal);
        end
    end
end