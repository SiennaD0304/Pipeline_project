# used wget (link to SRA file)
# afterwards I used fasterq-dump./ # with # being each SRA number 
# I did this for each SRA file and saved it to the directory folder generated with copying to GitHub 

#### Step 2 ####
# rule all for spades, finding the longest contig blast results 
# is using the SRA fastq numbers for the rule too 
sample = ["SRR5660030","SRR5660033","SRR5660044","SRR5660045"]
rule all:
    input:
        expand("results/spades_assemblies/{sample}", sample=sample),
        expand("results/longest_contigs/{sample}_contigs.txt", sample = sample),
        expand("results/blast_results/blast_{sample}.csv", sample = sample)

#getting genome from NCBI 
rule downloading_genome:  
    output:
        directory("dataset/GCF_000845245.1"), 

# the second line "datasets..." was provided by NCBI to get the genome information for Human herpesvirus 5
# adding cds and protein sequences to different output file to make it easier to find later 
    shell:
        """
        mkdir -p dataset/GCF_000845245.1
        datasets download genome accession GCF_000845245.1 --include cds,protein,genome --filename "dataset/ncbi_dataset.zip" 
        unzip -o "dataset/ncbi_dataset.zip" -d dataset/
        cp dataset/ncbi_dataset/data/GCF_000845245.1/cds_from_genomic.fna {output}
        cp dataset/ncbi_dataset/data/GCF_000845245.1/GCF_000845245.1_ViralProj14559_genomic.fna {output}
        cp dataset/ncbi_dataset/data/GCF_000845245.1/protein.faa {output}
    
        """
# make a list of all coding sequences within the genome and identify how many coding sequences there are 
rule make_list: 
    input:
        "dataset/GCF_000845245.1/cds_from_genomic.fna"
    output: 
        "dataset/output.fasta",
        "dataset/count_cds.txt"

# run command will let you run python inside snakefile! 
    run:
        #establishing variables and downloading biopython 
        from Bio import SeqIO
        seq = []
        id_list = []
        x = []
        seq_pair= []
        protein_list = []

        # making list of sequence names and sequence reads 
        for record in SeqIO.parse(input[0], "fasta"):
            id_list.append(record.id)
            seq.append(str(record.seq))

        #print(seq)
        for i in id_list: 
            x = list(i) # breaking up every letter into a list to find Y
            #print(x)
            for y in x:
                if y == "Y":
                    pos = x.index(y)
                    protein_name = (x[pos:]) # if protein name present get full name 
                    protein_name = ''.join(protein_name)
                    protein_name = "".join([">", protein_name]) #addin carrot 
                    protein_list.append(protein_name) #keeping protein names in a list
        count = 0
        for w in range(len(seq)): # getting 1st protein name paired with 1st dna seq ... 2nd name to 2nd seq
            seq_pairs = "".join([protein_list[w] +"\n", seq[w]])
            seq_pair.append(seq_pairs) #pairing each together
            count += 1 # adding to count for each pair made 

        #print(count)

        #outputting to sequence file for kallisto
        with open (output[0], "w") as outfile: 
            outfile.write("\n".join(seq_pair))

        with open (output[1], "w") as outfile:
            outfile.write(str(f"The HCMV genome (GCF_000845245.1) has {count} CDS."))

# making kallisto index 
rule kallisto:
    input: 
        "dataset/output.fasta" # the cds file 
    output:
        "dataset/kallisto_output.idx", #where I want index to be store 
    shell:
       "kallisto index -i {output} {input}" #command from slides 

#### STEP 3 ####
rule kallisto_quantify:
    input: 
        idx = "dataset/kallisto_output.idx", #index made from above rule 
        d1_2in_f = "SRR5660030_1.fastq", #forward and reverse SRA reads that we unzipped in the beginning 
        d1_2in_r = "SRR5660030_2.fastq",
        d1_6in_f = "SRR5660033_1.fastq",
        d1_6in_r = "SRR5660033_2.fastq",
        d2_2in_f = "SRR5660044_1.fastq",
        d2_2in_r = "SRR5660044_2.fastq",
        d2_6in_f = "SRR5660045_1.fastq",
        d2_6in_r = "SRR5660045_2.fastq"

    output:
        d1_2 = directory("results/SRR5660030"), #have one output folder within the results for each timepoint (ex. donor 1 day 2... )
        d1_6 = directory("results/SRR5660033"),
        d2_2 = directory("results/SRR5660044"),
        d2_6 = directory("results/SRR5660045")

    shell:
        # running kallisto: allows you to identify how many transcripts (from samples) you have for each section of a genome agains a refrence genome (index)
        # results in abundance reads for each sample 
        """
        kallisto quant -i {input.idx} -o {output.d1_2} -b 30 {input.d1_2in_f} {input.d1_2in_r} 
        kallisto quant -i {input.idx} -o {output.d1_6} -b 30 {input.d1_6in_f} {input.d1_6in_r}
        kallisto quant -i {input.idx} -o {output.d2_2} -b 30 {input.d2_2in_f} {input.d2_2in_r}
        kallisto quant -i {input.idx} -o {output.d2_6} -b 30 {input.d2_6in_f} {input.d2_6in_r}
        """
rule make_sleuth_table: 
    input: 
        d1_2 = "results/SRR5660030/abundance.h5",
        d1_6 = "results/SRR5660033/abundance.h5",
        d2_2 = "results/SRR5660044/abundance.h5",
        d2_6 = "results/SRR5660045/abundance.h5"
    output: 
       table = "results/sleuth.tsv"
    run: 
        #will use for sleuth in the next step 
        import pandas as pd #used pandas to make a table based off the abundance of trasncripts from kallisto 
        
        #making matrix format: "column_name" : info for each column 
        matrix_frame = {
        "sample":["SRR5660030","SRR5660033","SRR5660044", "SRR5660045"],
        "condition":["2dpi", "6dpi", "2dpi", "6dpi"],
        "path":[input.d1_2, input.d1_6, input.d2_2, input.d2_6]
        }

        sleuth_matrix = pd.DataFrame(matrix_frame)

        with open (output.table, "w") as outfile:
            outfile.write(str(sleuth_matrix))

#running sleuth with transcript abundance reads 
rule sleuth_run: 
    input: 
        table = "results/sleuth.tsv"
    output: 
        sleuth_out = "results/sleuth_output.txt"
    script:
        "sleuth.R" #sleuth script is attached within Github 


#### step 4 ####
#using the genome downloaded from rule 2 in snakefile 
rule bowtie2_genome_index: 
    input: 
        genome = "dataset/ncbi_dataset/data/GCF_000845245.1/GCF_000845245.1_ViralProj14559_genomic.fna",
    
    output: 
        directory("results/bowtie_data_hcmv")
    shell: 
        #write the directory before hand to make sure output is written beforehand 
       """
       mkdir -p {output} 
       bowtie2-build {input.genome} {output}/GCF_index
       """ 
       # I needed to do the /GCF_index to make sure that the files are written as GCF_index.(bt...) - it didn't work without this 

rule running_bowtie_2:
    #using the genome index made above with fastq forward and reverse reads downloaded in the beginning 
    input:
        genome_results = "results/bowtie_data_hcmv/GCF_index.1.bt2",
        d1_2f = "SRR5660030_1.fastq",
        d1_2r = "SRR5660030_2.fastq",
        d1_6f = "SRR5660033_1.fastq",
        d1_6r = "SRR5660033_2.fastq",
        d2_2f = "SRR5660044_1.fastq",
        d2_2r = "SRR5660044_2.fastq",
        d2_6f = "SRR5660045_1.fastq",
        d2_6r = "SRR5660045_2.fastq",

    output: 
        out_d1_2 = "results/bowtie_data_out/SRR5660030.sam",
        out_d1_6 = "results/bowtie_data_out/SRR5660033.sam",
        out_d2_2 = "results/bowtie_data_out/SRR5660044.sam",
        out_d2_6 = "results/bowtie_data_out/SRR5660045.sam"

        
    shell: 
        #used format from slides about bowtie two and the no-unal 
        #no-unal will get rid of any unpaired fastq files 
        """
        mkdir -p results/bowtie_data_out
        bowtie2 --quiet -x results/bowtie_data_hcmv/GCF_index -1 {input.d1_2f} -2 {input.d1_2r} -S {output.out_d1_2} --no-unal
        bowtie2 --quiet -x results/bowtie_data_hcmv/GCF_index -1 {input.d1_6f} -2 {input.d1_6r} -S {output.out_d1_6} --no-unal
        bowtie2 --quiet -x results/bowtie_data_hcmv/GCF_index -1 {input.d2_2f} -2 {input.d2_2r} -S {output.out_d2_2} --no-unal
        bowtie2 --quiet -x results/bowtie_data_hcmv/GCF_index -1 {input.d2_6f} -2 {input.d2_6r} -S {output.out_d2_6} --no-unal
        """

#converting the sam output of bowtie2 to a fastq 
rule sam_to_fastq:
    input: 
        "results/bowtie_data_out/SRR5660030.sam",
        "results/bowtie_data_out/SRR5660033.sam",
        "results/bowtie_data_out/SRR5660044.sam",
        "results/bowtie_data_out/SRR5660045.sam"
    output: 
        directory("results/sam_fastq_out")
    shell: 
        # made the directory first to have a file for snakefile to refrence back to 
        # used a for loop to change the name of each file 
        # base=$(basename $name .sam) I got from AI I couldn't figure out how to do it a different way 
        """
        mkdir -p {output}
        for name in {input}; do 
            base=$(basename $name .sam)
            samtools fastq $name -o {output}/$base.fastq
        done
        """
        #using samtools I go from stack overflow 

rule counting_reads:
    input:
        "SRR5660030_1.fastq",
        "SRR5660033_1.fastq",
        "SRR5660044_1.fastq",
        "SRR5660045_1.fastq",
        "results/sam_fastq_out/SRR5660030.fastq",
        "results/sam_fastq_out/SRR5660033.fastq",
        "results/sam_fastq_out/SRR5660044.fastq",
        "results/sam_fastq_out/SRR5660045.fastq"
    output:
        "results/bowtie_counts.txt"

    run:
        from Bio import SeqIO

        names = ["SRR5660030","SRR5660033","SRR5660044","SRR5660045" ]
        with open (output[0], "w") as outfile:

            for i in range(4):  #iterate start from 0-3 * bc non inclusive 
                #establishing counters ** inside loop so they reset ** 
                count_seq_1 = 0
                count_seq_2 = 0

                #starting at 0 for counter (aka before fastqs)
                for record in SeqIO.parse(input[i], "fastq"):
                    count_seq_1 += 1 

                #starting at 3 and iterating through counts (aka after fastq)
                for record in SeqIO.parse(input[i+4], "fastq"):
                    count_seq_2 += 1 

                #for single end reads (read_sum_before) I doubled it to account for the paired reads in bowtie 
                read_sum_before = count_seq_1*2 
                read_sum_after = count_seq_2 
                #having each round of input write to an outfile 
                outfile.write(str(f"Sample {names[i]} had {read_sum_before} read pairs before and {read_sum_after} read pairs after.\n"))

#### STEP 5 #### 
#used a rule all here because spades only wanted 1 input an output at a time and not 4 at once 
#goal was to generate contigs 
rule running_spades: 
    input:
        "results/sam_fastq_out/{sample}.fastq",

    output:
        directory("results/spades_assemblies/{sample}") #this line was used for the rule all so it would recognize the sample 

    shell:
        "spades.py -k 127 -t 2 --only-assembler -s {input} -o {output}" # code structure was from class, not sure if -t rule was acutally needed 
        # kmers were specificed in instructions 

# sorting contigs with python to find the longest contig for blast 
# used rule all again 
rule find_longest_contig: 
    input:
        "results/spades_assemblies/{sample}/contigs.fasta"

    output:
        "results/longest_contigs/{sample}_contigs.txt" #added this command to the rule all to use same sample list 

    run:
        #establishing variable 
        seq_len = []
        from Bio import SeqIO #importing biopython package (SeqIO)

        with open (output[0], "w") as outfile: # opened file first

                for record in SeqIO.parse(input[0], "fasta"):
                    seq_len.append(record.seq)
        
                longest_contig = sorted(seq_len, key = len, reverse = True) #sorting list with longest contig 1st 
                seq_len = [] # emptying vairable so it can be used for next sequence 

                outfile.write(str(longest_contig[0] + "\n")) #making 4 different files with the longest contig for each sequence 

#### STEP 6 ####
#making a local database from betaherpesvirinae genomes in NCBI 
rule blast_local_db:
    output:
        new_directory = directory("genome_local_db/blast_datasets"),
        tracing_directory = "genome_local_db/betaherpesvirinae.nsq"
    shell: 
        #the datasets line of code is from the framework provided by NCBI 
        #unzip tells where to open the file to 
        #makeblast framework is from slides
        """
        datasets download virus genome taxon betaherpesvirinae --include genome --filename ncbi_dataset.zip
        unzip -o ncbi_dataset.zip -d {output.new_directory}
        makeblastdb -in genome_local_db/blast_datasets/ncbi_dataset/data/genomic.fna -out genome_local_db/betaherpesvirinae -title betaherpesvirinae -dbtype nucl
        """

#running blast with rule all too 
rule do_blast: 
    input: 
        betaherpesvirinae_genome = "genome_local_db/betaherpesvirinae.nsq",
        sample_contig = "results/longest_contigs/{sample}_contigs.txt" # added this line into the rule all with extend as well 
    output:
        "results/blast_results/blast_{sample}.csv"
    shell: 
        """
        blastn -query {input.sample_contig} -db genome_local_db/betaherpesvirinae -outfmt "6 sacc pident length qstart qend sstart send bitscore evalue stitle" -max_target_seqs 5 -max_hsps 1 -out {output}
        """
        #had to generalize genome_local_db/betaherpesvirinae vs genome_local_db/betaherpesvirinae.nsq to make sure all files from folder are read 
        # max_target_seqs will only output the top 5 sequences

rule make_blast_outfile:
    # using a for loop again to make sure that the blast files are combined into one file for the PipelineReport
    output: 
        "results/blast_output.txt"

    shell:
        """
        list="SRR5660030 SRR5660033 SRR5660044 SRR5660045" 

        for blast_seq in $list; do 

            echo "$blast_seq" >> {output}
            cat results/blast_results/blast_$blast_seq.csv >> {output}
            echo >> {output}
        done
        sed  -i -e '1i sacc   pident   length   qstart   qend   sstart   send   bitscore   evalue   stitle' {output}
       """ # sed allows header to be added to the very first line of a file 
       #echo will add text to the file >> will append to output file 
       #cat will copy a file that exists to the output file
       #allows for all 4 blast output files to be read in as 1 

rule make_final_file: 
    # adding all of the output files into PipelineReport.txt 
    input: 
       cds_count_file = "dataset/count_cds.txt",
       sleuth_output = "results/sleuth_output.txt",
       bowtie_filtering_output = "results/bowtie_counts.txt",
       blast_output = "results/blast_output.txt"

    output:
        "PipelineReport.txt"

    shell:
        """
            cat {input.cds_count_file} >> {output}
            echo >> {output}
            echo >> {output}
            cat {input.sleuth_output} >> {output}
            echo >> {output}
            cat {input.bowtie_filtering_output} >> {output}
            echo >> {output}
            cat {input.blast_output} >> {output}
            echo >> {output}
        """        
        # echo >> output is there to add lines between each file 