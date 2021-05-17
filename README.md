# Metagenomics-lab
# **STRIM STRatification In Metagenomics**

## Using STRIM

### Cite Us
This code was developed as part of a project carried out during the course of Microbial Metagenomics (Molecular Biology Master Degree) at the University of Padova. The project was supervised by Prof. Stefano Campanaro and Dr. Arianna Basile.

## **Description**
The STRIM code is useful to calculate the abundances weighted for the occurrencies for each metabolic function (e.g. a function identified in the KEGG database)and it allows to define the stratification of the taxonomies for every different function. In other words, the stratification of the taxonomies means that the code enables the user to determine which is the contribution of each taxa at different taxonomic level for each specific function. 
The code takes as input:
* the abundances of each species in a community of a particular sample, calculated with checkM 
* the metabolic functions of each species of the community, annotated independently with eggNOG
* the taxonomic assignments containing the NCBI classification 
and it enables to interbreed these information, giving as output the stratification of the taxonomies contributing to each function.

## Requirements
In order to use the software, you will need:
* pandas: it is a library which allows you to perform data analysis and manipulation, particularly it is useful for manipulating numerical tables and time series.
* numpy: it is a library useful to manipulate large arrays and matrices and it offers high-level mathematical functions to work on these arrays.
* matplotlib.pyplot: it is a plotting library which enables to create static, animated or interactive plots and figures.
* os module: it enables the interaction between the user and the operating system.
* glob module: it is used to recover files or pathnames matching a specified pattern.
* fnmatch: supplies function to match files having a specific pattern or filter files having a specific pattern.
* collections: gives alternatives to create container data types such as list, tuple and dictionary.
* itertools: supplies functions that operate on iterators to produce complex iterators. An iterator is an object that contains a countable number of values that you can pass through, meaning that the vaues can be iterated upon. 

## Generation of the inputs and inputs required
The inputs needed will be:
* a file containing the abundances calculated with checkM
* a file containing the eggNOG annotations
* a file containing the taxonomic assignments with the NCBI classification
Pay attention: if you used different tools to calculate the abundances or the annotations, you have to modify the code with the correct inputs (e.g. in the code, when you find '.annotations' you will need to change with the format you used to save your data)

## How to use 
In order to use the code, the user will need to follow the subsequent steps:
* download the file strim.py
* open the file strim.py with a text editor and modify the paths, meaning the directories in which the files you are going to use are saved (e.g. every time you will have to save an output that is going to be generated 'examplefilename.to_csv', you will have to specify the path in which the file needs to be saved)
* check to have all the files required 
* change the columns names of the taxonomic assignment with the ones you have in your file containing the taxonomic assignments (e.g. column_list = ['average_Re1','average_Re2','average_Re3', 'average_Re4'] put the correct names between '')
* open the terminal and enter in the directory where the file .py is saved
* write 'python3 strim.py' and press the Enter key
* the program will ask you a series of parameters as a input:
  - the choice between KEGG orthologues and KEGG modules
  - the taxonomic level desired for the stratifcation step
* the program will elaborate the data and will return a message on the screen every time it will end the abundances calculation of a defined sample
* when the abundances calculation step is endend, the stratification step begins
* at the end of the strtatification step, the program will return a message on the screen, inviting the user to check the results obtained, opening the directory in which the files have been saved.

## Output
In the folder in which the user decided to save te outputs, the following output files should be found:
* two tabular files:
  - one defining the taxonomic assignment
  - one defining the abundances weighted fot the occurences for each metabolic function
  -  
a questo punto l'utente dovrà entrare nella cartella in cui ha salvato i dati e qui potrà trovare i due file tabulari, uno che definisce il taxonomic assignment e uno che è il file delle abbondanze, una serie di foto ognuna delle quali corrisponderà alle tot weighted abundances per ogni campioni. la stratificazione restituisce file tabulari uno per ogni kegg code preso in esame
