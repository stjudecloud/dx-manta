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
    mkdir inputs
    mv in/normal_bam/**/* inputs/
    if [[ ! -z "$normal_bam_index" ]]; then
      mv in/normal_bam_index/**/* inputs/
    fi

    CONFIG_CMD="configManta.py"

    # Check input BAMs
    if [[ -z "$tumor_bam" ]]; then
      echo "[*] No tumor sample. Assuming normal bams only"

      for bam in inputs/*.bam
      do
        echo "Checking normal bam file $bam"
        bam_name=`basename $bam`
        if [[ -z "inputs/$bam_name.bai" ]]; then
          echo "[*] No index file found for bam. Generating index..."
          samtools index $bam
        fi
        CONFIG_CMD="${CONFIG_CMD} --bam ${bam}"
      done
    else
      echo "[*] Tumor BAM detected. Will pair with one normal BAM sample."

      bams=(inputs/*.bam)
      echo "DEBUG: ${bams}"

      if [ "${bams[@]}" -gt 1]; then
        echo "[*] Found multiple normal bams, only one will be picked."
      fi

      bam=${bams[0]}
      echo "DEBUG: ${bam}"
      bam_name=`basename $bam`
      if [[ -z "inputs/$bam_name.bai" ]]; then
        echo "[*] No index file found for normal bam. Generating index..."
        samtools index $bam
      fi
      CONFIG_CMD="${CONFIG_CMD} --normalBam ${bam}"

      echo "[*] Checking tumor bam file $tumor_bam_name"
      mv in/tumor_bam/* inputs/
      mv in/tumor_bam_index/* inputs/
      if [[ -z "inputs/$tumor_bam_name.bai" ]]; then
        echo "[*] No tumor bam index file found. Generating index..."
        samtools index inputs/$tumor_bam_name
      fi

      CONFIG_CMD="${CONFIG_CMD} --tumorBam inputs/${$tumor_bam_name}"
    fi

    # Setup Reference Fasta
    if [[ -z "$reference_fasta" ]]; then
      echo "[*] No reference FASTA found. Will choose a default hg38 reference."
      if samtools view -H inputs/$normal_bam_name | grep EBV; then
        # EBV found
        echo "[*] Using GRCh38_no_alt"
        dx download project-F5444K89PZxXjBqVJ3Pp79B4:file-F89b7z09gk6zYbY98vJvVX37 -o inputs/
        dx download project-F5444K89PZxXjBqVJ3Pp79B4:file-F89b7z09gk6YQPfJ93kb9K54 -o inputs/
        CONFIG_CMD="${CONFIG_CMD} --referenceFasta inputs/GRCh38_no_alt.fa"
      else
        # No EBV found
        echo "[*] Using hg38m1x"
        dx download project-F5444K89PZxXjBqVJ3Pp79B4:file-Fy4BkZj9PZxyvJQb586JXgZP -o inputs/
        dx download project-F5444K89PZxXjBqVJ3Pp79B4:file-Fy4Bp1Q9PZxjgVZ1KqV0KG1Z -o inputs/
        CONFIG_CMD="${CONFIG_CMD} --referenceFasta inputs/hg38m1x.fa"
      fi
    else
      echo "[*] Reference FASTA found."
      mv $reference_fasta_path inputs/
      if [[ -z "$reference_fasta_index" ]]; then
	echo "[*] No reference FASTA index found. Indexing now..."
	samtools faidx inputs/$reference_fasta_name
      else
	mv $reference_fasta_index_path inputs/
      fi
      CONFIG_CMD="${CONFIG_CMD} --referenceFasta inputs/$reference_fasta_name"
    fi

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
    mv run/results/stats out/statistics_outputs
    dx-upload-all-outputs
}

