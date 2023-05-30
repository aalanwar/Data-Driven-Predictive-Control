## Robust Data-Driven Predictive Control using Reachability Analysis
<br /> 
<br /> 
This repo cotains the code for our two papers:<br /> 
[1] Amr Alanwar*, Yvonne Stürz*, Karl Johansson "Robust Data-Driven Predictive Control using Reachability Analysis" European Journal of Control <br />
<br /> 
<br /> 
We present a robust data-driven control scheme for unknown linear systems with a bounded process and measurement noise. Instead of depending on a system model as in traditional predictive control, a controller utilizing a data-driven reachable region is proposed. The data-driven reachable regions are based on a matrix zonotope recursion and are computed based on only noisy input-output data of the system's trajectory. We assume measurement and process noise which are contained in bounded sets. While we assume knowledge of these bounds, no knowledge about the statistical properties of the noise is assumed. Our proposed scheme guarantees robust constraint satisfaction under measurement and process noise, which is essential in safety-critical applications.<br />

<br /> <br />
<p align="center">
<img
src="Figures/reachmpc.png"
raw=true
alt="Subject Pronouns"
width=500
/>
</p>
[2] Mahsa Farjadnia*, Amr Alanwar, Muhammad Umar B. Niazi, Marco Molinari, and Karl Henrik Johansson. "Robust Data-Driven Predictive Control of Unknown Nonlinear Systems using Reachability Analysis." Accepted to European Control Conference (ECC) 2023. <br />
<br /> 
<br /> 
We present a robust data-driven predictive control approach for unknown nonlinear systems in the presence of bounded process and measurement noise. By using the past noisy input-output data in the learning phase, we propose a novel method to over-approximate reachable sets of an unknown nonlinear system. Then, we propose a data-driven predictive control approach to compute safe and robust control policies from noisy online data. The constraints are guaranteed in the control phase with robust safety margins through the effective use of the predicted output reachable set obtained in the learning phase.   <br />
<br /> 
<br />
<br /> 

## Running 
### paper [1] <br />
<br />
1- Download [MPT](https://www.mpt3.org) and install [mosek](https://www.mosek.com/products/academic-licenses/) toolboxs.<br />
2- Add MPT folder and subfolders to the Matlab path.  <br />
3- Add the whole folder of this repo and subfolders to the Matlab path.  <br />
4- run ZPC.m.<br />
5- run Robust_MPC_polytopes.m.<br />
6- run plotPolyZono.m <br />

### paper [2] <br />
<br />
1- Download and install the [MPT](https://www.mpt3.org) and [mosek](https://www.mosek.com/products/academic-licenses/) toolboxes. <br />
2- Add the MPT folder and its subfolders to the Matlab path. <br />
3- Include the entire folder of this repository, along with its subfolders, in the Matlab path. <br />
4- Execute the NZPC_Predictive_Control.m file. <br />
<br />
<br />

## Ack
Note that portion of this code is from CORA toolbox and from Felix Gruber, and Matthias Althoff "Scalable Robust Model Predictive Control for Linear Sampled-Data Systems"
We acknowledge the efforts by Paul George (UWaterloo) in fixing one Bug in the code. Many thanks!


Our paper Bibtex is as follows:<br />
```
@article{alanwar2022robust,
  title={Robust data-driven predictive control using reachability analysis},
  author={Alanwar, Amr and St{\"u}rz, Yvonne and Johansson, Karl Henrik},
  journal={European Journal of Control},
  pages={100666},
  year={2022},
  publisher={Elsevier}
}
```
```
@article{farjadnia2022robust,
  title={Robust Data-Driven Predictive Control of Unknown Nonlinear Systems using Reachability Analysis},
  author={Farjadnia, Mahsa and Alanwar, Amr and Niazi, Muhammad Umar B and Molinari, Marco and Johansson, Karl Henrik},
  journal={arXiv preprint arXiv:2211.05867, Accepted at European Control Conference (ECC)},
  year={2023}
}
```
