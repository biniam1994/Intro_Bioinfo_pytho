# Intro_Bioinfo_python

### **FastQC_Parser**

The purpose of the application is to parse fastqc file as input and generate appropriate graphs for each fastqc module and txt file to specified directory. The application is fully implimented using python programming language. The application will run from command line terminals in unix/linux, mac and windows command prompt. 

# Running the pragramm

## **Step 1 - Installation**
## **Dependencies**
Install this dependencies.
 * Python 3.8.3 
 * pip 20.2.4
 * And the following python packages or libraries:
 
 | Package | Version | 
 |-------|-------|
 | `Seaborn` | 0.4.6 |   
 | `Matplotlib` | 2.1 |  
 | `pandas` | 6.0.86 | 
 | `numpy` |  
 | `argparse` | 0.21 | 
 | `os` | 0.21 | 
 | `errno` | 0.21 | 
## **Step 2 - Running python programm from command line**
run below command to see for help with -h flag and know the usage of argument to be passed to obtaine specific result
```
$ python Assigment.py -h
```
Example we can pass analyze the basic stastics by specifing -m or -module = module name (see help to check short notation of modules),
-i input fastqc file, -o output directory or folder.
```
$ python Assigment.py -i fastqc_data1.txt -o Result -m BS
```
The same here in this command as well we are only changing the -m module to `ALL`, which will run all the module and output the file to 
appropraite directory
```
$ python Assigment.py -i fastqc_data1.txt -o ALL_results -m ALL
```

```
$ python Assigment.py -i fastqc_data1.txt -o Filter_report -m Filter
```
