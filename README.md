# EC552-HW-1
## Description
Genetic circuit design program for modifying custom Eco1C1G1T1.UCF.json and Eco1C1G1T1.input.json files

## Instructions
#### 1. Install python 3.8+ and pip
#### 2. Install Docker
`docker pull cidarlab/cello-dnacompiler:latest`
#### 3. Install poetry
`poetry install` 

`poetry run python main.py`
#### 4. Cello file structure

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
├── poetry.lock
├── pyproject.toml
└── README.md
```

