# QC Workflow
- `wgs_qc_wf.wdl`: WDL QC workflow for WGS data

> [DONE] test fastp
```bash
dx run /apps/wdl_wf/wgs_qc_wf/fastp_wf \
    -istage-common.fastq_r1=file-GX0Gv3002k8vFY5xGb9469xj \
    -istage-common.fastq_r2=file-GX0Gv3802k8vYyKj90gJx397 \
    -istage-common.prefix="UDN921066-P" \
    -y --brief --folder /gcarvalho_test/qc_wf/test_fastp \
    --name UDN921066-P_fastp
```

> [DONE] test picard
```bash
dx run /apps/wdl_wf/wgs_qc_wf/picard_qc \
    -istage-common.bam_file=file-GX0QyVj02k8z7zbJJKkpqk63 \
    -istage-common.bam_index=file-GX0QyZ802k8qQq5BBYPFxp7B \
    -istage-common.fasta=project-BQpp3Y804Y0xbyG4GJPQ01xv:file-FF2vqv007JZyg5vFFBYb0gJZ \
    -istage-common.fasta_dict=project-BQpp3Y804Y0xbyG4GJPQ01xv:file-GFz5xf00Bqx2j79G4q4F5jXV \
    -istage-common.fasta_fai=project-BQpp3Y804Y0xbyG4GJPQ01xv:file-FFJx1P80XJyP87xzF632jqqQ \
    -istage-common.prefix="UDN921066-P_chr1" \
    -y --brief --folder /gcarvalho_test/qc_wf/test_picard \
    --name UDN921066-P_chr1_picard_qc
```

> [TODO] test qualimap
```bash
dx run /apps/wdl_wf/wgs_qc_wf/qualimap_bamqc_wf \
    -istage-common.bam_file=file-GX0QyVj02k8z7zbJJKkpqk63 \
    -istage-common.bam_index=file-GX0QyZ802k8qQq5BBYPFxp7B \
    -istage-common.prefix="UDN921066-P" \
    -y --brief --folder /gcarvalho_test/qc_wf/test_qualimap \
    --name UDN921066-P_qualimap
```

> [DONE] test multiqc
```bash
dx run /apps/wdl_wf/wgs_qc_wf/multiqc_wf \
    -istats_files=file-GX17xj80jgJvZxJzxJF71QzY \
    -istats_files=file-GX17xj80jgJvYyKj90gKbX6z \
    -istats_files=file-GX183780JQG6jJG9YVBk4xbG \
    -istats_files=file-GX1837Q0JQGBq9169KP5F0Xq \
    -istats_files=file-GX183780JQGP4J5j47X5Yx0z \
    -istats_files=file-GX183780JQG8F1yvkb15zK1k \
    -istats_files=file-GX183780JQG7qBf4G6fJbp5x \
    -istats_files=file-GX183780JQG6P0GPFkQPXzf2 \
    -istats_files=file-GX183780JQGJFQkbZbXpYPQK \
    -istats_files=file-GX183780JQG4B13B1x15Y6Bb \
    -istats_files=file-GX1837Q0JQGFQq5BBYPGKZzG \
    -istats_files=file-GX183780JQG06Xx034JXGkKp \
    -istats_files=file-GX1897j0zzfBqgbXB7Qk0vzF \
    -iprefix="UDN921066-P_filtered" \
    -y --brief --folder /gcarvalho_test/qc_wf/test_multiqc \
    --name UDN921066-P_filtered_multiqc
```

> test everything
```bash
dx run /apps/wdl_wf/wgs_qc_wf/qc_wf \
    -istage-common.bam_file=file-GX0QyVj02k8z7zbJJKkpqk63 \
    -istage-common.bam_index=file-GX0QyZ802k8qQq5BBYPFxp7B \
    -istage-common.fasta=project-BQpp3Y804Y0xbyG4GJPQ01xv:file-FF2vqv007JZyg5vFFBYb0gJZ \
    -istage-common.fasta_dict=project-BQpp3Y804Y0xbyG4GJPQ01xv:file-GFz5xf00Bqx2j79G4q4F5jXV \
    -istage-common.fasta_fai=project-BQpp3Y804Y0xbyG4GJPQ01xv:file-FFJx1P80XJyP87xzF632jqqQ \
    -istage-common.fastq_r1=file-GX0Gv3002k8vFY5xGb9469xj \
    -istage-common.fastq_r2=file-GX0Gv3802k8vYyKj90gJx397 \
    -istage-common.prefix="UDN921066-P_filtered" \
    -y --brief --folder /gcarvalho_test/qc_wf/test_qc \
    --name UDN921066-P_filtered_qc
```

> test multi sample for multiqc
```shell
dx run /apps/wdl_wf/wgs_qc_wf/find_data_multiqc_wf \
    -istage-common.dir_path=/Analysis/hg38_udn/UDN921066/filtered_bam/ \
    -istage-common.prefix="UDN921066-P_UDN085245-F_UDN501807-M" \
    -y --brief --folder /Analysis/hg38_udn/UDN921066/filtered_bam/ \
    --name UDN921066-P_UDN085245-F_UDN501807-M_multiqc
```
---
