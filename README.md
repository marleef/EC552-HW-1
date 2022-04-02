# EC552-HW-1
### Authors
Drew Gross and Marlee Feltham

## Description
Genetic circuit design program for modifying custom Eco1C1G1T1.UCF.json and Eco1C1G1T1.input.json files

## Instructions
#### 1. Install python 3.10 and pip
#### 2. Install Docker
`docker pull cidarlab/cello-dnacompiler:latest`
<!-- #### 3. Install poetry
`poetry install` 

`poetry run python main.py` -->
#### 3. Cello file structure

`git clone https://github.com/CIDARLAB/homework1-template`

```
├── input
│   ├── and.v
│   ├── Eco1C1G1T1.input.json
│   ├── Eco1C1G1T1.output.json
│   ├── Eco1C1G1T1.UCF.json
│   ├── nand.v
│   ├── options.csv
│   ├── struct.v
│   └── xor.v
├── main.py
├── output
└── README.md
```

## Suggested Reading
1. Vaidyanathan, P., Der, B. S., Bhatia, S., Roehner, N., Silva, R., Voigt, C. A., & Densmore, D. (2015). A Framework for Genetic Logic Synthesis. Proceedings of the IEEE,
103(11), 2196-2207.
(This Paper should help understand the basic framework, terms and definitions of Logic Synthesis applied for Genetic circuits. Specifically, Sections
1,2,3 as well as Fig 2 will be helpful)
2. Yaman, F., Bhatia, S., Adler, A., Densmore, D., & Beal, J. (2012). Automated selection of synthetic biology parts for genetic regulatory networks. ACS synthetic biology,
1(8), 332-344.
(The Results section as well as Fig 5 will be useful)
3. Stanton, B. C., Nielsen, A. A., Tamsir, A., Clancy, K., Peterson, T., & Voigt, C. A.
(2014). Genomic mining of prokaryotic repressors for orthogonal logic gates. Nature
chemical biology, 10(2), 99-105.
