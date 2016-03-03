# ** Optimisation of Asset Management Maintenance  **

##Aim
The aim of the project is to classify crossarms into either wood or steel from the power pole images given from United Energy’s image data set. By doing so, this will enhance the asset management maintenance by reducing time and cost for people to manually look through them. Various computer vision and image processing techniques were explored to fully understand the approaches to be taken and steps to ultimately achieve the desired goal. R programming language was used. 

##System Overview

![System Overview](https://github.com/UnitedEnergy/Automated_Pole_Inspector/blob/master/System_Flow.png)

The above flowchart depicts the system overview of this project. 
Image pre-processing technique was used by using Component Labelling. This script and documentations can be found in the Component Labelling folder of this project. The purpose of this is to identify whether there is a crossarm present to reduce the amount of images that should be processed or looked at. The reason for this is there are so many images in the data set, with a lot of ‘useless’ images in which there are no crossarms present at all. 
The second technique is the Feature Based Approach where it can predict the crossarm in the image as steel or wood. More of this approach can be found in the Feature Based Approach folder. 
There are also two PowerPoint files, one which is a presentation, the other a poster. This can be looked at to get a quick overall understanding of the project and the techniques that were incorporated to succeeding the goal. 

