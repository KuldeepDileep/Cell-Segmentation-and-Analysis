main();

%Main function:
function main()
    input_img = imread('polymersomes.tif');
    processed_img = preprocess(input_img);  %This is task 1 where preprocessing of the input image is performed.
    img_edge = detectEdge(processed_img);   %This detects edges.
    clear_img = clear(img_edge);            %This clears non-cellular components from image (method 1).
    seg_img = segment(input_img);           %This is another method to perform segmentation of cell. 
    cell_analysis(clear_img);%This extracts properties and performs regional property analysis
    cell_analysis(seg_img); %This extracts properties and performs regional property analysis
end

%Task1:
function output = preprocess(input)
    figure
    img = im2double(input);
    imshow(img(:,:,3))
    title('Original Image');
    
    %Noise Filtering:
    figure
    filter = fspecial('gaussian', [8 8]);   %Applying 8x8 gaussian filter
    filter_img = imfilter(img, filter);
    imshow(filter_img(:,:,3))
    title('Gaussian Filtered Image');

    figure
    filter2 = fspecial('average', [8 8]);   %Applying 8x8 averaging filter
    filter_img2 = imfilter(img(:,:,3), filter2, 'replicate');
    imshow(filter_img2)
    title('Averaging filtered Image');

    figure
    img2 = img(:,:,3);
    patch = imcrop(img2,[170, 35, 50 50]);
    patchVar = std2(patch)^2;
    DoS = 2*patchVar;
    filter_img3 = imbilatfilt(img2, DoS);   %Applying an edge-preserving Gaussian bilateral filter where DoS is the degree of smoothing
    subplot 121
    imshow(patch)
    title('Patch')
    subplot 122
    imshow(filter_img3)
    title('Gaussian bilateral filtered Image')

    figure
    x = wiener2(img(:,:,1), [8 8]);
    y = wiener2(img(:,:,2), [8 8]);
    z = wiener2(img(:,:,3), [8 8]);
    Inew(:,:,1) = x;                        %Applying a pixel-wise 2-D 8x8 adaptive low-pass Wiener filter
    Inew(:,:,2) = y;
    Inew(:,:,3) = z;
    imshow(Inew)
    title('2-D adaptive noise-removal filtering')

    %Sharpening:
    figure
    sharp_img = imsharpen(filter_img2,'Radius',2,'Amount',2);   %Sharpening the averaging filtered image
    imshow(sharp_img)
    title('Sharpened Image (Averaged Image)')

    figure
    output = imsharpen(Inew(:,:,2),'Radius',2,'Amount',2);      %Sharpening the 2-D apative filtered image
    imshow(output)
    title('Sharpened Image (Wiener)')
end

%Task 2:
function output = detectEdge(input)
    %Edge Detection:
    figure
    output = edge(input, 'canny', 0.28, 4);         %Applying canny edge detection on the sharpened image
    imshow(output)
    title('Edge Image')
end

%Task 3:
function output = clear(input)
    figure
    filled_img = imfill(input,'holes');             %Filling the regions in the edge image
    imshow(filled_img);
    title('Filled Image');

    figure;
    output = bwareaopen(filled_img, 300, 4);        %This removes the objects with number of pixels less than 300 to give us the final image with 6/7 cells segmented this is also known as area opening operation
    imshow(output);
    title('Final Segmented Image - 1');
end

%Task 4:
function output = segment(input)
    figure
    img = im2double(input);
    imshow(img(:,:,3))
    title('Original Image');
    %Step 1:
    figure
    CLAHE_img = adapthisteq(img(:,:,3));            %This performs adaptive histogram equalization
    imshow(CLAHE_img)
    title('Contrast-limited adaptive histogram equalization (CLAHE)')
    %Step 2:
    figure;
    struct_element = strel('disk', 5);
    marker = imerode(CLAHE_img, struct_element);    %This erodes(shrinks) the regions in the equalized image
    imshow(marker)
    title('Eroded Image')
    %Step 3:
    figure;
    reconstruct_img = imreconstruct(marker, CLAHE_img); %This reconstructs the eroded image
    reconstruct_marker = imadjust(reconstruct_img);     
    imshow(reconstruct_marker)
    title('Reconstruct Marker')
    %Step 4:
    figure
    stack(:, :, 1) = reconstruct_marker;
    %Step 5:
    sigma = [1 1.5 2 2.5 3 3.5];
    for i = 6                                           %This stackes up more 6 more channels filtered at different sigmas
        filter = fspecial('gaussian',[83 83],sigma(1,i));
        stack(:, :, i+1) = imfilter(reconstruct_marker, filter, 'replicate');   %Filtered channels being added to the stack
    end
    %Step 6:
    double_stack = im2double(stack);
    %Step 7: 
    ones_stack = ones(size(double_stack), 'double');
    %Step 8:
    new_stack = double_stack + ones_stack;
    %Step 9:
    std_dev = std(new_stack, 0, 3);                     %Returns standard deviation of the new stack  
    %Step 10:
    std_dev = std_dev./mean(new_stack, 3);              
    std_dev = double(std_dev);
    %Step 11: 
    std_dev = std_dev-min(std_dev(:));
    %Step 12:
    std_dev = 255-im2uint8(std_dev);    
    %Step 13:
    std_dev = im2double(std_dev);
    std_dev = im2double(reconstruct_marker)+std_dev;
    %Step 14:
    std_dev = mat2gray(std_dev);
    thresh = 0.85*graythresh(std_dev);
    binary_img = imcomplement(imbinarize(std_dev,thresh));  %We complement the resulting image to make background black
    imshow(binary_img);

    figure;
    filled_img = imfill(binary_img,'holes');
    imshow(filled_img);
    title('Filled Image');

    figure;
    labelled_Image = bwlabel(filled_img);
    stats = regionprops(labelled_Image,'area','Perimeter'); %Gives area and perimeters of all objects

    areas = [stats.Area];                           %Gives us circular regions in the image
    perimeters = [stats.Perimeter];
    boundary = perimeters.^2 ./ (4*pi*areas);
    regions = find(boundary <1.9);
    binary_Image = ismember(labelled_Image, regions) > 0;
    imshow(binary_Image);
    title('Circular Regions')

    figure                                          %This performs opening that is erosion and the dilation to only preserve circular bodies (size of structuring element is set such that is completely eats up non-circular regions 
    struct_element3 = strel('disk', 10);
    open_img = imopen(binary_Image,struct_element3);
    imshow(open_img)
    title('Opening')

    figure
    struct_element4 = strel('disk', 12);    
    open_img2 = imopen(binary_Image,struct_element4);   %This performs opening that is erosion and the dilation to further process in the image by only giving the 7 desired cell bodies and removing the extra at bottom right  
    imshow(open_img2)
    title('Opening')

    figure;
    size_of_window = 15;                                %This smoothens the boundaries
    kernel = ones(size_of_window) / size_of_window ^ 2;
    blurry_image = conv2(single(open_img2), kernel, 'same');
    output = blurry_image > 0.5; 
    imshow(output);
    title('Final Segmented Image - 2');
end

%Task 5:
function cell_analysis(input)
    img = imread('polymersomes.tif');
    label_image = bwlabel(input);                           %Labels all the bodies in the cell segmented image
    cell_array = [];                                        %To hold all the properties
    reg_prop = regionprops(label_image,'area','Perimeter'); %Gives regional property such as area and perimeter of every region

    for i=1:length(unique(label_image))                     %This is executed for the number of objects in the image
        [x, y] = find(label_image == i);                    %We find coordinates where pixel value = i 
        mask = zeros(size(input));                          

        for j = 1:length([x,y])
            mask(x(j),y(j)) = 1;                            %This creates a mask with only one body/region at a time
        end

        img_new = im2double(img(:, :, 1)).*mask;            %We apply it on the original image to acquire the desired region
        reg_prop2 = regionprops(mask);
        binary_img = imbinarize(img_new,0.8);

        figure
        imshow(binary_img)
        title('Dots Segmented')
        
        binary_image_label = logical(binary_img);   
        reg_prop3 = regionprops(binary_image_label,'area');

        num_dots = size(reg_prop3);
        if num_dots>0                                       %Regions with dots are visited to acquire properties such as number of dots in that region and ratio
            cell_array = [cell_array, [num_dots(1), sum([reg_prop3.Area],2)/reg_prop(i).Area]];
        end
    end

    num = length(cell_array)/2
    fprintf('Number of cells in the original image = %s \n', int2str(num));

    num = [];
    ratios = [];
    for i=1:length(cell_array)                  %This loop generates two separate arrays for number of dots in each region and ratio of area of each dot to total area of cell
        if rem(i,2)==0
            ratios = [ratios, cell_array(i)];
        else
            num = [num, cell_array(i)];
        end
    end

    fprintf('Number of dots in each cell = %s \n', int2str(num));
    fprintf('Ratio between the area of dot and cell:');
    disp(ratios)
end
