# Effloc: An efficient locating algorithm for mass-occurrence biological patterns with FM-Index
~version 0.1 (<u>20240320</u>)~   
version 0.2 (<u>20241223</u>)   
ToDO: optimize code structure, add friendly comments 

# What is it?
Effloc is an efficient locating algorithm for mass-occurrence biological patterns with FM-Index. 

# How to use it?
Effloc consists of three components, index building, pattern generating, and pattern locating. You should first build the FM-Index with the reference genome, and then perform the locating process.   
![0e8370ed8ffe521689206ec159e99e12](https://github.com/Lilu-guo/Effloc/assets/23703069/8d599869-5c0f-4894-82ec-98d10b5303de)


## Step I. Install
  1. Download (or clone) the source code form https://github.com/Lilu-guo/Effloc
  2. Compile the source code.
     ```shell
     cd ./Effloc   
     make
     ```
## Step II. Build FM-index
  1. Run the shell command:
     ```shell
     ./effloc <Follow the prompts to enter the genome file name>
     ```
## Step III. Generate pattern
  1. Run the shell command:
     ```shell
     ./effloc <Follow the prompts to enter the pattern set name>
     ```
## Step IV. Locate pattern
  1. Run the shell command:
     ```shell
     ./effloc <Follow the prompts to enter the pattern set name and index file name>
     ```
![68eae00ab1a97bead284329cc4a1ba3d](https://github.com/Lilu-guo/Effloc/assets/23703069/7a1ee058-47a3-4f10-936d-772b46aa2f58)

# Feedback
Please report bugs to Email: guolilu@stu.xidian.edu.cn or jgxygll@163.com if any questions or suggestions.   
Your feedback and test results are welcome.

# License   
Effloc is available under the MIT license.   
