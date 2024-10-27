# VFA T1 Mapping Using Derivative-Free Optimization Techniques

This project was completed for the course MTH 574, Spring'24, and experiments with alternative derivative-free optimization methods for VFA T1 mapping in MRI, inspired by the NOVIFAST method.

## Overview
This repository contains MATLAB implementations for estimating T1 maps in magnetic resonance imaging (MRI) using various derivative free optimization techniques. It includes functions for implicit filtering, Nelder-Mead optimization, and a fast optimization algorithm known as Novifast. The scripts provide a framework to add noise to the synthetic signals and assess the impact of different signal-to-noise ratios (SNRs) on the estimated T1 maps.

## File Descriptions

- **demo_implicit_filtering_image.m**: Demonstrates the implicit filtering optimization method for estimating T1 maps. 
- **demo_ndelder_mead_image.m**: Demonstrates the use of the Nelder-Mead optimization method for estimating T1 maps.
- **demo_novifast_image.m**: Provides a demonstration of the Novifast optimization technique for estimating T1 maps.
- **demo_novifast_SNR.m**: Illustrates how to assess the performance of the Novifast method under various SNR conditions.
- **if_image.m**: Implements the implicit filtering algorithm for T1 estimation.
- **implicit_filtering_optimization.m**: Contains the optimization logic for the implicit filtering method.
- **model_based_optimization.m**: Provides a framework for model-based optimization approaches.
- **nelder_mead_image.m**: Implements the Nelder-Mead optimization algorithm for T1 estimation.
- **nelder_mead_optimization.m**: Contains the optimization logic for the Nelder-Mead algorithm.
- **nelderVsnovifast.m**: Compares the performance of the Nelder-Mead and Novifast optimization methods.

## Requirements
- MATLAB (version R2016b or later recommended)
- Image Processing Toolbox (for any image processing functions)

## Usage Instructions

1. **Load Data**: Before running any of the demo scripts, ensure you have the necessary data files, such as `slicePhantom.mat`, in the appropriate directory (./data).

2. **Running Demo Scripts**:
   - To run a specific demo script, open MATLAB and navigate to the directory containing the scripts.
   - Run the script by typing its name in the command window. For example:
     ```matlab
     demo_novifast_image
     ```
   - Each demo script will generate visual outputs and may save results in specified formats.

3. **Modifying Parameters**:
   - You can modify key parameters within the demo scripts, such as SNR levels and optimization options, to experiment with different settings.

4. **Reviewing Results**: 
   - Each demo script will generate figures and error metrics that shows the effectiveness and accuracy of the optimization methods in estimating T1 maps.

## Example
To run the Novifast demo with a specific configuration:
```matlab
% Navigate to the folder containing the scripts
cd 'path/to/your/scripts'

% Run the Novifast demo
demo_novifast_image
```

## References
## References

1. Maurovich-Horvat, P., Ferencik, M., Voros, S., Merkely, B., & Hoffmann, U. (2018). Comprehensive plaque assessment by coronary CT angiography. *Nature Reviews Cardiology*, 15(6), 303-313. Retrieved from [https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6277233/](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6277233/)
   
2. Wilm, B. J., Barmet, C., Pruessmann, K. P., & Kasper, L. (2018). Novifast: A Fast Algorithm for Accurate and Precise VFA MRI. Retrieved from [https://www.mathworks.com/matlabcentral/fileexchange/67815](https://www.mathworks.com/matlabcentral/fileexchange/67815)

3. Author(s). (Year). Title of the article. *Journal Name*, Volume(Issue), Page numbers. Retrieved from [https://pubmed.ncbi.nlm.nih.gov/22807160](https://pubmed.ncbi.nlm.nih.gov/22807160)

4. Author(s). (Year). Title of the article. Retrieved from [https://www.math.uci.edu/~qnie/Publications/NumericalOptimization.pdf](https://www.math.uci.edu/~qnie/Publications/NumericalOptimization.pdf)
