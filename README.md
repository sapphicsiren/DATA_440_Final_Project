# DATA_440_Final_Project
DATA 440 Final Project

## **Background**
### *Objective*
A program that can be utilized to classify and label theoretically generated conformational isomers for various ProXxx dipeptides. The rationale behind this project is that the manual labeling of these isomers is something already done in my research lab, and it is a process that could definitely be semi-automated. 

For this first stage of the project, I am focused on just classifying the dipeptides as either cis or trans isomers. 

### *Project Context*
In the directory titles 'slides' there is a slideshow pdf which includes some background chemistry information along with an explanation for how this project is set-up. 

## **The Data**

The raw data will not be made publically available as it contains all of the theoretical research performed by my lab. Instead, a cleaned version of the data is uploaded, and it includes the cartesian coordinate data for each atom and the classification label of the isomer (c/t). 

## **Current Functionality**

Currently my program is able to read, preprocess, visualize, and classify each isomer as cis or trans.  

To run the code:
- Clone the repository
- Install plotly.express if you do not have it using
  !pip install plotly.express
- Open function_playground.ipynb and run blocks

### *Playground Breakdown*
Important Parameters:
- dipep: this is the variable containing the name of the ProXxx dipeptide that is being analyzed. Within the function_playground, the specific dipeptide being looked at is ProAla.
- dipep_conformer: this is the variable containing the name of the specific conformer of the ProXxx dipeptide we are looking at. Within the first portion of the function_playground, this is proala0054, a conformer of ProAla.

Important Functions: Note that each of these functions contains numerous smaller functions that are defined in more detail in the .py files
- vaa.make_conf_tables(dipeptide, dipeptide_conformer): from the 'visualization_and_algorithm_visible.py' file, this function is used to read the coordinate information for the specified conformer and output the data as a sql table.
- vaa.conf_visualization(dipeptide_conformer): from the 'visualization_and_algorithm_visible.py' file, this function is used to visualize the specified conformer in a 3-dimensional scatterplot. Note that the function returns a subsetted dataframe of the coordinates excluding all hydrogen atoms and a figure that visualzies the conformer without hydrogen atoms.
- vaa.c_t_isomer(dipeptide_conformer): from the 'visualization_and_algorithm_visible.py' file, this function is used to classify the specified conformer as either a cis or a trans isomer. The algorithm does so by comparing the angle between the two bonds of interest.
- vaa.validate_ct_labeling(dipeptide, dipeptide_conformer): from the 'visualization_and_algorithm_visible.py' file, this function is used to validate the classification algorithm by comparing the predicted label to the actual label.

Results: As shown in the final block in the notebook, the classification algorithm has 100% accuracy when it comes to labeling all of the ProAla conformers as either cis or trans. 
