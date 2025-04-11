# SAM (Sequence Alignment/Map) Processing

Sam.py is a Python script designed to process SAM (Sequence Alignment/Map) files, commonly used in bioinformatics to store large nucleotide sequence alignments. This tool facilitates the extraction and analysis of alignment data, providing users with a streamlined approach to handle SAM files.

## Features

SAM File Processing: Parses SAM files to extract relevant alignment information.

Data Analysis: Computes statistics such as alignment counts, mapping quality distributions, and more.

Output Generation: Produces summary reports or filtered SAM files based on user-defined criteria.

## Installation

Clone the Repository:
```
git clone https://github.com/BowerH/Sam.py.git
cd Sam.py
```

## Usage

The script can be executed via the command line:

```
python sam.py -i <input.sam> -o <output.txt> [options]
```

## Command-Line Arguments

-i, --input: Path to the input SAM file. (Required)

-o, --output: Path to the output file. (Required)

-q, --min-quality: Minimum mapping quality to include in the analysis. (Optional)

-f, --flag: SAM flag to filter alignments. (Optional)

#### Example

```
python sam.py -i eg_data/sample.sam -o results/summary.txt -q 30
```

This command processes sample.sam, filters alignments with a mapping quality of at least 30, and writes the summary to summary.txt.

## Output

The script generates an output file containing:

Total number of alignments.

Number of alignments passing filters.

Distribution of mapping qualities.

Other relevant statistics.












