# Autosvp
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
   
### Initializing Meta-Space and Running Lumpy-sv with Multi-Threading- Python version 3.6 or higher
 ```bash
   bash sample_histoandmeastd.sh
   python mulit_2lumpy.py  meanstdev.csv 
   ```
### Gradient-Free Parameter Optimization Using Gaussian Processes and Bayesian Optimization
 ```bash
   
   python GP.py  your.csv 
   ```
### To generate metadat：
   your.csv + namefeature.csv
###  Training a meta-model, in the context of machine learning, involves creating a model that can learn from the outputs or the performance of other models.then get multi_target_regression_model.pth
 ```bash
   
   python AUTOsvp.py 
   ```
### Testing model with recommended parameters 
```bash
   
   python model.py \Autosvp\samplecopy\1.bam
   print
    prediction_dict = {
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
### Testing model in LUMPY-SV
 ```bash
   
   python test.py yourpath/PreorRecallorf1.csv
   ```

###  Predicting neoantigens by .vcf in manual & Auto  parameters recommended framework 
 THE computational method termed NeoSV, which incorporates SV annotation, protein fragmentation, and MHC binding prediction together, to predict SV-derived neoantigens. https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-03005-9
  ```bash
   
    neosv -vf test.sv.vcf -hf test.hla.txt -np /path/to/netmhcpan -o test -p test -r 75
  ```
### Evaluating the accuracy of a VCF (Variant Call Format) file in manual & Auto  parameters recommended framework 
  ```bash
    python eva.py result.csv
  ```
   
