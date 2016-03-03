# **Automated Pole Inspector**

##Aim
The aim of the project is to classify crossarms into either wood or steel from the power pole images given from United Energy’s image data set. By doing so, this will enhance the asset management maintenance by reducing time and cost for people to manually look through them. Various computer vision and image processing techniques were explored to fully understand the approaches to be taken and steps to ultimately achieve the desired goal. R programming language was used. 

##System Overview

![System Overview](https://github.com/UnitedEnergy/Automated_Pole_Inspector/blob/master/System_Flow.png)

The above flowchart depicts the system overview of this project. 

Image pre-processing technique was used by using Component Labelling. This script and documentations can be found in the Component Labelling folder of this project. The purpose of this is to identify whether there is a crossarm present to reduce the amount of images that should be processed or looked at. The reason for this is there are so many images in the data set, with a lot of ‘useless’ images in which there are no crossarms present at all. 

The second technique is the Feature Based Approach where it can predict the crossarm in the image as steel or wood. More of this approach can be found in the Feature Based Approach folder. 

There are also two PowerPoint files, one which is a presentation, the other a poster. This can be looked at to get a quick overall understanding of the project and the techniques that were incorporated to succeeding the goal. 

##Component Labelling

###Overview
The aim of this procedure is to determine whether an image contain a crossarm or not. In order to achieve this, we employed different techniques that will transform the image, identify the components that are present on it and select those ones that tend to look like crossarms sections just to give a final picture’s classification pending on the filtered regions analysis.

The full procedure is managed by one script that takes the pictures from a given working directory and processes them to obtain a list of those images and its segmentation that is stored in an output file (csv). 

###Code Features and Phases

####Directories and main elements
The script expects a directory that contains the images that are about to be analyzed (working directory) and it allows to add different storage paths for the results obtained on different stages of the procedure. The most relevant points on this section are:

-Directories definition.

-Size delineation for reduced images. The variable “j” usually was tested with value of 500 (500 x 500 pixels) but it can be adjusted.

-Variable “output” that will store the picture classification is created.


####Image Pre-processing
In this section the following elements are implemented:

-List of all files (pictures) present on the working directory.

-Reduced version (“img”) of the original image.

-Edge detection is performed using the gradient (dilate – erode) on the blue layer of the resized picture. The blue layer was picked as a grayscale version of the original file due to its higher contrast levels between some color shades that tend to be problematic to handle in this image set such as green areas, sky, among others.
 
####Component Labelling
This module contains the main processes of the program, from creating the elements required to apply this technique and performing the two-pass algorithm till obtaining the areas to be described and selected in the next modules. The main sections of the process are described as it follows:

####Preparations
-A label matrix (“labels”) is created; this has the same size as the reduced image and will storage the area’s label that every pixel has been associated to.

-An equivalence table (“areas”) is defined in order to record the existent relations between areas that are part of the same picture’s element but have been labelled on different regions.

-The threshold value (“v”) that identifies if the image sections are related. This value is adjustable pending on the sensibility wanted to perceive changes in the pixel values. We decided to use 0.0225 since when applying the gradient to identify the edges on the image it is relative easy to notice the changes across the different elements present on the image, this value represents a change of 2.25% in the whole scale of possible values and it’s as well the 50% of the average standard deviation obtained in successfully labelled pictures.

####Two-pass (based) algorithm 
-The image is scanned (from left top to right bottom corner) and the pixels are labelled according on its relation with the ones that are located at the left and up of the focused pixel.

-The areas are analyzed and any related regions are identified even if they got different labels on the previous step.

-The labels obtained are exported to a csv file (optional) as a backup for analysis purposes.

####Area Description 
-An area catalogue (“cat”) is defined to storage size and scope of every region created.

-The minimum and maximum number of pixels that an area should contain in order to be considered as potential crossarm section are defined. These values (“minval” and “maxval”) are adjustable and are dependent on the size of the reduced image (“j”). The actual value consider that usually a crossarm section does not cover more than 20% of the image but has always at least 0.2% of the total pixels enclosed by the resized picture.

-The values to store in the catalogue are calculated and stored on “cat”.

####Area Selection and Filtering
-Evaluates if the regions founded match the size required to be considered as potential crossarm sections discarding all those that are outside that stablished range of values.

-Once the areas that have a reasonable size are kept from dismissing, they are analyzed to verify if the shape, density and location are consistent with the features present on the crossarm pictures.

-After discarding those sections that didn’t match the previous set of statements, the remaining areas are being filtered on a copy of the resized image in order to be able to recognize the original colors on them.

-For removing more irrelevant areas, color thresholding is applied to this filtered version of the image so we can recognize certain color shades that are distant to the metallic and wood crossarm color values.

####Component Labelling (2nd time)
After finishing the preceding filtering processes, component labelling is performed again so the remaining pixels can be reassigned according to the original edges of the image. This second version of the technique is almost equal to the one already performed to the image, containing the following phases:

-Preparations.

-Two-pass (based) algorithm. For this occasion, the adjustable thresholding value (“v”) was redefined to .02 since the previous color filtering processes lessened the original sections edges and a more sensible scan worked better for more accurate regions recreation.

-Area description.

-Area Selection and Filtering. In this module, some of the selection statements have been tailored since it is expected to have “smaller” areas than the ones obtained on the first labelling, so the region filtering considers this fact to still be consistent with the proportional filtering applied throughout the different sections of this program. In addition to this, some of the color filters have been discarded since there is no point on look for pixels with values already removed from the processing image.

####Picture Classification
This unit creates a copy of the filtered area labels, so all the pixels enclosed on these regions have the same value (the area label’s number) having all the discarded areas as background with pixel values equals to 0. The method implemented follows the next steps to complete the image segmentation using this regions analysis:

-The regions are picked one at a time and it is verified that the selected region is not labelled already as background in order to proceed to the analysis. It is proven too that the picked area matches the required minimal length (“j/20”) defined to be considered as potential crossarm section. This value is adjustable too.

-Once the area matched the first two validations, the matrix “column” is created. This matrix will store the number of pixels on each column that are associated to the selected area, in order to count how many columns match the required minimum height (given the number of pixels on a column) for that given area. That minimum height is given by “j/33”, which is adjustable.

-The rows enclosed on the selected region are scanned and the pixels on each column associated to the region picked are counted and stored on “column”.

-The columns matching that minimum height value (“j/33”) are counted to verify if this number matches the minimum of columns required to consider that area a crossarm section. This number of required columns is defined by “j/10” (adjustable) and alongside “j/33” was defined to seek the minimum values needed to preserve a size relevant area with the minimum length and height to consider it as potential section to classify the whole image.

-If the picked section does not match these validations, the area is redefined as background and the loop will continue until the first area that matches those criteria appears (meaning there is a crossarm on the picture) or when there is no remaining areas to check.

-After this loop is completed, the picture is classified into the variable “tag” and stored alongside the name of the image (“filename”) into the matrix “output”, which will accumulate the list of images processed with their resulting classification.

-Finally, when all images have been classified, the content of the “output” variable is exported into a csv file to preserve the results obtained.

##Feature Based Approach


###Code Overview – How to use
Functions.R, totalfunc.R, and folderfunc.R must be run first. They contains the relevant libraries and some functions that are needed by the other files. Note that some of the libraries may need to be installed (using install.packages) the first time the code is run.

Put some photos in the training directory. Images for training need to be copied into the \\calculation files\\training folder.

Training.R contains code for the first part of training. The first three lines contain variables that can/should be edited:

-Overall.direc is the directory containing the photos to be processed. It must also contain the ‘calculation files’ folder (see Directory structure section)

-No_cores indicates the number of available cores the code is permitted to use

-No_train is the number of photos to be used for training the neural network. This number must be less than the number of photos present in the training directory – the remaining photos will be used for testing and verification.

Also of importance is line 29. By default commented out, it generates a csv file where the user can manually label training photos. This must be done when new training photos are added, and it is important that any ‘type names’ do not contain spaces. This step does not need to be repeated on subsequent runs and hence can be commented out.

Trainingpart2.R can be run after Training.R and as long as a .csv file (called ‘type labels’) with labels for training photos exists. This part of the code trains the neural network, and outputs a percentage accuracy as calculated on the test set. This section can be run multiple times if a higher accuracy is desired, however warnings may be generated on subsequent runs (these can safely be ignored)

The workspace should be saved at this point so that the trained neural network can be recovered later if required.

Testing.R then uses the neural network from the previous code to classify unknown images. The outputs of this code appear as .csv files in the ‘outputs’ folder.

####Required Directory Structure
Overall.direc refers to the folder containing the pole photos (these need to be grouped in folders labelled by pole number). It also assumes that there is more than one photo per pole – if this is not the case it may give an error, which can be easily overcome by placing a dummy photo in the relevant folder.

The ‘calculation files’ .zip folder also needs to be unzipped and copied into this same directory. It contains files that the code needs to run. This is also where the ‘type labels.csv’ file appears, which can be edited by the user to manually label training photos.

![File Structure](https://github.com/UnitedEnergy/Automated_Pole_Inspector/blob/master/File_Structure.png)

####Functions.R
Functions.R contains a variety of user-defined functions for use by the other files.

#####Padimage(I, p)

Takes an image (I) and inserts a buffer of p pixels around the edges, such that filters/fourier transforms can be applied to the image without distorting the edges.

#####Blur_pic(pic, sigma, ksize)

Applies a Gaussian blur of strength sigma and filter width ksize to an input picture pic.

#####Pixsample(x,y)
Performs bilinear interpolation on the current working image at given x and y values. This is used for downsampling images.

#####Patch(sigma, pic, x_cor, y_cor)
Takes a picture (pic) and returns a patch of that image, centered around x_cor and y_cor, created with a gaussian window of strength sigma.

#####X_edge(pic)
Takes a picture (pic) and uses a sobel filter to extract the x-direction edges.

#####Y_edge(pic)
Takes a picture (pic) and uses a sobel filter to extract the y-direction edges.

#####Gradients(pic)
Takes an input picture (pic) and returns both the magnitude and orientation of the gradients present in the image.

#####Folderfunc and totalfunc
Refer to training and testing sections, respectively.

####Training.R
Training.R begins by defining relevant directories, and then sets up totalfunc() to process training images across multiple cores. Totalfunc works to create a csv file containing descriptors of the keypoints in a training image. This is done by the following basic steps:

1. Sequential Gaussian blurs

2. Difference of gaussians (highlighting key points)

3. Repeat 1-2 over several scales

4. Non-maximal suppression (discard keypoints if the neighbour is more significant)

5. Take patch around each keypoint

6. Rotate to standard orientation

7. Calculate gradient-based descriptor

8. Write to output file

Once descriptors have been calculated for all training images, these descriptors are then grouped into 250 pre-calculated clusters. This clustering is done in a two-layer tree structure for efficiency, and each point is only calculated as a ‘good match’ for a given cluster if it is significantly closer to that centre than to any other, helping to discard poor matches.

The overall output is then a single .csv file containing each image represented as a 250 element vector of the number of points it contained belonging to each cluster.

####Trainingpart2.R

Training part 2 takes the cluster output file from training part 1 and combines it with user-entered labels from the ‘type labels’ files to train a neural network. Although training part 1 created an output for every image in the training directory, training part 2 uses only a subset of these (of size no_train) to train the neural network – the remaining images are kept separate and used only for testing.

Training part 2 hence begins by separating the data into two – a test set and a training set. It then imports the ‘correct’ user-defined labels for each image, and uses the training data to create a neural network. The neural network currently has a single hidden layer of 102 nodes, however this can be modified for optimisation.

Once trained, the network is tested on the ‘test’ images, and a percentage accuracy returned based on these.

At this point, the workspace should be saved so that the trained neural network can be retrieved again later.

####Testing.R

Testing.R follows much the same process as the training. Folderfunc is used to process multiple photos at a time – each core is used to process an entire pole at once (a folder’s worth of photos). Descriptors are created for each photo, these descriptors are clustered, and then processed through the neural network to output predictions of type of photo.

The output .csv files are written to the outputs directory and contain are structured as follows:

Pole number | Photo name | Steel.crossarm | Useless | Wood.crossarm
------------|------------|----------------|---------|--------------
Pole 1|Photo1.jpg|0.000215|2.76E-05|0.98
Pole 1|Photo2.jpg|0.004194|0.973|0.0673
Pole 1|Photo3.jpg|0.001304|0.09923|5.35E-06

The readings in the latter three columns give an indication of how likely it is that the photo belongs to that category – with 1 being the maximum.

####Combineouts.R

Combineouts can be used once testing is completed. It merely combines all the separate output files (one for each pole) into one .csv file called ‘combined outputs’, which can be found in the ‘calculation files’ folder.
