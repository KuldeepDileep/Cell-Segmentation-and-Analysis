# Cell Image Segmentation and Analysis

Biological cells generally exhibit heterogeneity requiring them to be analyzed on a single cell basis. In order to do that, individual cells need to be extracted from the image. This is a complicated task as the cells may overlap or touch each other and the segmentation not only means separating cells from background but also from each other. Once the individual cells are extracted then any desirable process could be applied on them to perform biological analysis. Here we take up an example scenario where we want to segment cellular as well as sub-cellular regions from image and then apply some random processing and analysis on them.

## Cell Image Segmentation:
### Task 1:
Take the image 'polymersomes.tif' and preprocess it, e.g., to convert it to appropriate data type, to enhance and remove any noise, and finally sharpen it.
### Task 2: 
Apply appropriate edge detection method to extract the edges of the cell body as is depicted in the Figure 1b. Once the edge detection is performed we can use it to segment the regions of interest, that is, the cell body. The edges can be connected if they are not already and flood fill operation (use help to study function for it) on the resulting image gives the segmented cells regions as shown in Figure 1c.
### Task 3: 
The resulting image has non-cellular components and false edges due to imperfect segmentation as visible in Figure 1c. Use a combination of morphological operations to perform filtering to get an image as shown in Figure 1d. Use help to study functions like bwareaopen, imclearborder for getting an idea if they could be applicable here. For example consider removing objects with pixel area less than say 200-400 pixels.
### Task 4: 
Apply cell image segmentation method described in the attached paper page 4. Compare the segmentation results from the two methods.

## Cell Image Analysis (Spots detection, counting and density):
### Task 5: 
Label (bwlabel) the cell regions for individually accessing them in performing single cell analysis. Using the value of label, create a row vector to contain number of spots on each cell. Also create cell array to hold regional properties ‘structure’ for each cell. Now, in a loop (using labels), access/extract the individual cell regions from the labeled image such that the image contains only that region (use Logical operator). Extract regional properties (regionprops) of the extracted region. Now, use extracted object image as mask to extract the corresponding region/object of the original intensity image. Threshold that image, setting a threshold, e.g., between 75%-90% of maximum value, will help extract further features, e.g., spots within the cells. Labeling the thresholded image (which is actually spots in cells) would give the number of spots in that particular cellular region. Calculate the size(or area)-number of spots ratio of the cells, e.g., to highlight their health.
### Task 6: 
Perform Task5 using the other segmentation result and comment on the accuracy of the results.

**Note: Read report uploaded with the code to visualize the results of each task.**
