# ACSE9_Independent-Research-Project

<!-- PROJECT SHIELDS -->

<!-- TABLE OF CONTENTS -->

<details open="open">
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
    </li>
    <li>
      <a href="#prerequisites">Prerequisites</a>
    </li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgements">Acknowledgements</a></li>
  </ol>
</details>

<!-- ABOUT THE PROJECT -->

## About The Project

The project contains ku used to compress and predict urban flow field models to effectively predict fluid dynamics simulations. 
The code was constructed for the master's thesis of Imperial College London.


<!-- GETTING STARTED -->

```
.
├── POD                  #POD coefficient method
│   ├── dataset          #Test data set (vtu file)
│   └── process          #.py file
├── Compress_AE          #Compression via AAE
│   ├── dataset          #
│   ├── model            #Saved model
│   └── notebook         #Notebook for training AAE
├── prediction           #Make predictions through the Adversarial network
│   ├── dataset          #
│   ├── model            #Saved model
│   └── notebook         #Notebook for training Adversarial network 
├── regular              #Apply on regular grid
│   ├── dataset          #regular data set
│   └── notebook         #Notebook used for calculation
├── result               #
│   ├── Animation        #Generated animation
│   └── VTU              #Generated VTU file 
├── preprocessing        # 
│   └── generate         # Preprocessing py files
└── tests                # Package tests
    └── data             # Test datasets(seen and unseen dataset)
```

Please note that for this project, the test is mainly done in the form of a global run in the entire build software in the colab notebook. These notebooks can be found in the above document. 
The data set path needs to be modified according to the actual situation.

## Prerequisites

* Python 2.7
* Tensorflow
* (Recommended) GPU with CUDA 11

## Reproduction of reported results

Since the compression and prediction models need to be trained for a long time, the trained model is also uploaded to the "model" folder. Because the data set is too large to upload to github, only a part of the vtu file is included as an example in the test folder. You can generate the data set for training through the .py file in the preprocessing folder. If you need a complete test vtu data, please contact me.

## Contact

* xaingqi.Liu chessliuxx@gmail.com

<!-- ACKNOWLEDGEMENTS -->
## Acknowledgements

* Prof. Christopher Pain
* Dr. Claire Heaney
* Royal School of Mines, Imperial College London
