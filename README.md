# Double Tagged Amplicon Demultiplexing

This software demultiplex experiments from an illumina run with double tagged amplicon. 
Attention: double tag don't mean double index.
The tags are present at the 5' border of the primers. They are not included in the illumina adapter.  
  
For more details please see:  
Microbiome Profiling by Illumina Sequencing of Combinatorial Sequence-Tagged PCR Products, Plos one, October 2010  
Gregory B. Gloor, Ruben Hummelen, Jean M. Macklaim, Russell J. Dickson, Andrew D. Fernandes, Roderick MacPhee, Gregor Reid  


## Usage

### Inputs

* experiments in a CSV file: The software need a CSV file containing all the experiments. This file must contain at least 4 fields names run, sample, forward and reverse. Run correspond to your illumina run id ; sample, the name of yout sample in the run ; forward and reverse, the names of your tagged primers.  
  
Here an example with 2 illumina runs on two sample with 2 primers and 3 tags:  

    run,sample,forward,reverse
    JAN_2017,st_1,PR1-A,PR2-B
    JAN_2017,st_2,PR1-B,PR2-C
    FEB_2017,st_1,PR1-A,PR2-B
    FEB_2017,st_2,PR1-B,PR2-C

* Primers FASTA file: To complete the CSV, a FASTA file with all the tagged primers is required. Each primer is represented by its name in header and its sequence.  
  
Here an example corresponding to the previous CSV example (tagsize: 4bp):  

    >PR1-A
    ACCTGCCTAGCGTYG
    >PR1-B
    GAATGCCTAGCGTYG
    >PR2-B
    GAATCTYCAAATCGG
    >PR2-C
    ACTACTYCAAATCGG

* R1/R2 FASTQ files: The paired-end FASTQ files outputed by the illumina sequencing machine.

### Outputs

For each experiment, the software will output 2 FASTQ files named run_sample_fwd.fastq and run_sample_rev.fastq. These files will contain all the paired forward and reverse reads.  
  
To not fill your current directory with many FASTQ files, you can define an output directory in the command line.

### Command line

    ./dtd -r1 <r1_filename.fastq> -r2 <r2_filename.fastq> -o <oligos_filename.fasta> -e <experiments.csv> [-d output_directory] [-m]
    
-m : Output unasigned reads in fastq files.

## Download and compile

The software is designed to compile on a unix system (linux or macOS).  

```bash
git clone git@github.com:yoann-dufresne/DoubleTagDemultiplexer.git DTD  
cd DTD/  
make  
```
