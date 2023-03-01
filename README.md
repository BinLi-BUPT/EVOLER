<a href="https://doi.org/10.5281/zenodo.7688653"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.7688653.svg" alt="DOI"></a>


This file provides an instruction for the operation software to reproduce all the figures/results of our work “Machine Learning-enabled
Globally Guaranteed Evolutionary Computation”.
1.	Operation software of evolutionary computing 
1)	Software:  MATLAB R2019a  
The source code has been tested on MATLAB R2019a and it could work normally. Note that, however, the other MATLAB versions may have unexpected problems, probably due to different specifications or versions of system functions. If the operational software (MATLAB R2019a) has been correctly installed, it requires no further action.  
 
2)	Operational system:    Windows 11 (64)
3)	Source Code Files:    
 
Note: There are 4 folders in the code repository, corresponding to the reproduction results of Figs. 2, 3, 4, 5 in the manuscript. 
Each folder contains the source code files (both our method and the comparative methods), the required data and the obtained results (data and figures). Please refer to each single Readme. doc file in every folder for the details information on source code. 
4)	Necessary Toolbox Installation: Tensor toolbox
  （1）：Unzip the tensor_toolbox folder under the Matlab folder paths
（2）：Enter the following code at the command line：
>> addpath('F:\matlab\tensor_toolbox')
>> addpath('F:\matlab\tensor_toolbox\met')  >> savepath 
>> doc sptensor
5)	Add the MATLAB Path to LUMERICAL:
(1)	Open the LUMERICAL; and then click the ‘Help->Matlab Integration status’
 
(2)	Add the MATLAB path in the following position with red box.
 
(3)	Open MATLAB as the administrator, and enter the following code at the command line：
>> !matlab /regserver
   
2.	Operation software of Finite Difference Time Domain (FDTD) simulation 
1)	Software:  LUMERICAL 2018a  
For the FDTD simulation of the nano-photonics inverse design, we use the LUMERICAL 2018a to obtain the resulting transmission response of specific nano-photonics structures. Note that, however, this software should be combined with MATLAB R2019a in the analysis. For example, the source code is running on MATLAB R2019a to obtain the potential solutions, while LUMERICAL 2018a evaluates the fitness of the provided solution. That is to say, they should be interacted in the whole optimization process.
Even so, in our EVOLER method, we have also made the best to minimize the complex interaction between MATLAB and LUMERICAL. For example, we have also provided the dataset derived via the LUMERICAL software, i.e., the high-dimensional tensor, which serves as the simple fitness function in the first stage of structured sampling, without running LUMERICAL frequently. 
  
2)	Operational system:    Windows 11 (64)
3)	Source Files:    
 
Note: There are 2 folders in the code sub-repository, corresponding to the reproduction results of Fig. 5 in the manuscript, i.e., in both the double transmittance band (Figs. 5d and 5e) and the single transmittance band (Figs. 5f and 5g). Please refer to the corresponding Read_me. doc file for the detailed information on this source code.
