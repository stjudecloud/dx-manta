#!/bin/bash

main() {
    set -exo pipefail

    # Install packages
    sudo apt-get update -qq
    sudo apt-get install -qq bzip2 gcc g++ make python zlib1g-dev samtools

    echo "Value of normal_bam: '$normal_bam'"
    echo "Value of normal_bam_index: '$normal_bam_index'"
    echo "Value of tumor_bam: '$tumor_bam'"
    echo "Value of tumor_bam_index: '$tumor_bam_index'"
    echo "Value of reference_fasta: '$reference_fasta'"
    echo "Value of reference_fasta_index: '$reference_fasta_index'"
    echo "Value of output_contig: '$output_contig'"

    echo "[*] Downloading all inputs"
    dx-download-all-inputs --parallel

    # Move input bams to inputs directory
    echo $HOME
    mkdir inputs
    mv $HOME/in/normal_bam/**/* inputs/
    mv $HOME/in/normal_bam_index/**/* inputs/

    CONFIG_CMD="configManta.py" # --bam input.bam --referenceFasta reference.fa --runDir run

    # Check inputs
    if [[ -z "$tumor_bam" ]]; then
      echo "[*] No tumor sample. Assuming normal bams only"
      for bam in $HOME/inputs/*.bam
      do
        echo "Checking normal bam file $bam"
        bam_name=`basename $bam`
        if [[ -z "$HOME/inputs/$bam_name.bai" ]]; then
          echo "[*] No index file found for bam. Generating index..."
          samtools index $bam
        fi
        CONFIG_CMD="${CONFIG_CMD} --bam ${bam}"
      done
    else
      echo "[*] Tumor BAM detected. Will pair with one normal BAM sample."

      normal_bams=($HOME/inputs/*.bam)
      echo "DEBUG: ${normal_bams}"

      if [ "${normal_bams[@]}" -gt 1]; then
        echo "[*] Found multiple normal bams, only one will be picked."
      fi

      normal_bam=${normal_bams[0]}
      echo "DEBUG: ${normal_bam}"
      normal_bam_name=`basename $bam`
      if [[ -z "$HOME/inputs/$normal_bam_name.bai" ]]; then
        echo "[*] No index file found for normal bam. Generating index..."
        samtools index $normal_bam
      fi
      CONFIG_CMD="${CONFIG_CMD} --normalBam ${normal_bam}"

      tumor_bam_name=`basename $HOME/in/tumor_bam/*.bam`
      echo "[*] Checking tumor bam file $tumor_bam_name"
      mv $HOME/in/tumor_bam/* inputs/
      mv $HOME/in/tumor_bam_index/* inputs/
      if [[ -z "inputs/$tumor_bam_name.bai" ]]; then
        echo "[*] No tumor bam index file found. Generating index..."
        samtools index $HOME/inputs/$tumor_bam_name
      fi

      CONFIG_CMD="${CONFIG_CMD} --tumorBam inputs/${$tumor_bam_name}"
    fi

    # Setup Reference Files
    mv $reference_fasta_path inputs/
    if [[ -z "$reference_fasta_index" ]]; then
      echo "[*] No reference FASTA index found. Indexing now..."
      samtools faidx inputs/$reference_fasta_name
    else
      mv $reference_fasta_index_path inputs/
    fi
    CONFIG_CMD="${CONFIG_CMD} --referenceFasta inputs/$reference_fasta_name"

    # Add runDir option
    mkdir run
    CONFIG_CMD="${CONFIG_CMD} --runDir run"

    # Add outputContig option
    if [ "$output_contig" = true ]; then
      CONFIG_CMD="${CONFIG_CMD} --outputContig"
    fi

    # Install manta
    echo "[*] Installing manta..."
    wget https://github.com/Illumina/manta/releases/download/v1.6.0/manta-1.6.0.centos6_x86_64.tar.bz2
    tar -xjf manta-1.6.0.centos6_x86_64.tar.bz2

    # Add manta to path
    export PATH=`pwd`/manta-1.6.0.centos6_x86_64/bin:$PATH

    echo "[*] Running Manta config..."
    echo $CONFIG_CMD
    # eval $CONFIG_CMD
    RESULT=$($CONFIG_CMD)
    echo $RESULT

    echo "[*] Running Manta workflow..."
    run/runWorkflow.py

    # Move and upload outputs
    mkdir out
    mv run/results/variants out/variant_outputs
    mv run/results/statistics out/statistics_outputs
    dx-upload-all-outputs

    exit 1
}

