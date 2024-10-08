# SVSwiftStabilizer
## ABSTRACT
Recent studies emphasized the vital roles of structural variations (SVs) in cancer diagnosis and treatments. Although many state-of-the-art variant callers have been developed during the past decade, their performance often exhibits significant fluctuation across heterogeneity tumor samples, which directly disturbs the biomarkers for treatments, e.g. introducing errors on identifying neoantigens from SVs. Actually, most variant callers are equipped multiple user-defined parameters to fit for heterogeneous samples; however, this function is rarely used due to two difficulties: computational complexity from adjusting interdependent parameters and the theorical instability for manually tuning strategies. To address these two challenges, we proposed SVSwiftStabilizer(3Ser), an automated parameter recommendation framework employing meta-learning with Bayesian optimization module. This framework utilizes Bayesian optimization to efficiently search for optimal parameters by learning from past evaluations, ensuring stability in detection performance across diverse samples. It also incorporates meta-learning to derive broad insights into parameter performance from extensive data, creating a powerful predictive meta-model. This strategy enables effective navigation to optimal or suboptimal solutions, significantly decreasing the time complexity associated with parameter interdependence in high-dimensional spaces, addressing the issue of NP-completeness. Benchmark results demonstrated that 3Ser delivers a stable and high-speed parameter strategy, notably lowering the false positive and false negative rates in the prediction of SV-derived neoantigens, as well as proving its effectiveness in the analysis of count-based biomarkers. This emphasizes its potential to advance the identification of cancer-specific neoantigens.

## Usage
### Simulated Data files Generation

Simulate fq data using GSDcreator by following the method described in the paper available on IEEE Xplore.

**Reference:** 
[Simulation method for fq data](https://ieeexplore.ieee.org/abstract/document/8983192)

**Requirements:**
- Python version 3.6 or higher

**Steps to Generate Simulated Data:**

1. **Generate Simulation Scripts**
   Use `A_stableCallerPaperSimFlow.py` to generate `sh` scripts for simulating 10,000 samples:

   ```bash
   python A_stableCallerPaperSimFlow.py -o sh_files
   ```
### Simulated Data Generation with GSDcreator

Generate simulation scripts for 10,000 samples, install GSDcreator, and run the simulations as follows:
**Requirements:**
- Python version 2.7
1. **Generate Simulation Scripts:**
   To create simulation scripts for 10,000 samples, use the following command:

2. **Install GSDcreator:**
Install GSDcreator and import the `sh_files`.

   ```bash
    dos2unix 10000run.py
    chmod +x 10000run.py
   ```
4. **Prepare the Simulation Script:**
Convert the script to Unix format and make it executable:

5. **Run the Simulation:**
Ensure Python version 2.7 then execute the script:
    ```bash
    python 10000run.py
    ```
6.**To distribute the converted `.vcf` files into each simulation folder：
Move your `shinvcf.py` script to the Python 2.7 environment's `bin` directory:**
   ```bash
     mv noinsertshinvcf.py /yourpath/py2env/bin/
     chmod +x /home/cloudam/.conda/envs/py2env/bin/shinvcf.py
     dos2unix /home/cloudam/.conda/envs/py2env/bin/shinvcf.py
     shinvcf.py A_stableCallerPaperSimFlowShell.sh base.vcf
```
7.To split the samples into segments ranging from thousands (kilobases) to millions (megabases) of base pairs, use the following approach:
- Python version 3.6 or higher.Splitting Samples with Multithreading Support:
   ```bash
   python split_bv32.py
   ```
### Extracting Sample Meta-Features
- Python version 3.6 or higher
   ```bash
   python feature.py-o namefeature.csv
   ```
   
### Initializing Meta-Space with Multi-Threading- Python version 3.6 or higher
 ```bash
   
   python /SVfolder/main.py
   ```
###  Training a meta-model, in the context of machine learning, involves creating a model that can learn from the outputs or the performance of other models.then get multi_target_regression_model.pth
 ```bash
   
   python AUTOsvp.py 
   ```
### Testing model with recommended parameters 
```bash
   
   python model.py \Autosvp\samplecopy\1.bam
   print
    prediction_dict = {
         "w"
        "msw": ,
        "tt": ,
        "back_distance": ,
        "min_mapping_threshold":],
        "min_clip":,
        "read_length":,
        "min_non_overlap":,
        "discordant_z": 
    }
   ```
### Testing DELLY, LUMPY, Manta, BreakDancer, Pindel, MetaSV, SvABA, and SVstabilizer
 ```bash
   
   python /SVfolder/vsworkflow/callerworkflow.py
   ```


   
