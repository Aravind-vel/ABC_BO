# Adaptive Boundary Constraint in Bayesian Optimization

**"Adaptive boundary constraint in Bayesian optimization: a general strategy to prevent futile experiments in complex reaction optimization"**  

Aravind Senthil Vel, Julian Spils, Daniel Cortés-Borda, François-Xavier Felpin  

📄 Published in *Reaction Chemistry & Engineering*, DOI: 10.1039/D5RE00161G 

🔗 [Access the paper here](https://doi.org/10.1039/D5RE00161G) 

This work introduces a strategy to avoid futile experiments in reaction optimization—i.e., experimental conditions that are theoretically found not to improve the current best objective value.  

This is achieved by incorporating a simple constraint into the optimization problem:

> **Max. achievable objective > current best objective found**

This constraint can be handled by Bayesian Optimization. They are known constraints, meaning the constraint violation can be validated before performing the reaction. Thus, the strategy can be used with any BO package that supports known constraints.

In this work, we used MATLAB's built-in function for Bayesian Optimization, bayesopt. 

As explained in the paper, we also perform boundary reduction for certain variables in each iteration to further avoid futile experiments. This step is optional; the constraint alone is sufficient to prevent futile experiments. However, reducing the boundaries may facilitate the algorithm more effectively identifying the remaining promising regions within the search space. 

---

## 📦 Repository Contents

**Optimization folder** contains the class files needed to define and solve an optimisation problem using either standard Bayesian Optimisation (BO) or ABC-BO, implemented with MATLAB’s built-in bayesopt function. The code is structured in an ask-tell manner, which are more suitable for reaction optimisation tasks. Since MATLAB only supports integer variables (for BO) by default, our code allows the user to specify variables as discrete numeric values (e.g. 1, 3, 8, 30 min). Constraints can also be defined over these discrete values.

**Examples** folder contains the scripts used to run both the in silico and real experimental case studies discussed in the paper.

**Results** folder contains the results of the case studies.  

---

## 🚀 How to Use

### 1. **Requirements**
- MATLAB 
- Optimization Toolbox
- Statistics and Machine Learning Toolbox

### 2. **Run an example**

For demo, refer to scripts in the **Examples** folder:  
`ISR_1_optimization.m`, `ISR_2_optimization.m`, `ISR_3_optimization.m`, `BHC_optimization.m` 

When running examples don't forget to add the **Optimization folder** to path:

