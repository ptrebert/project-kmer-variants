# project-kmer-variants

Conda environments defined under `environment/conda`, the Snakemake base environment is defined in the file `conda_smkenv.yml`.

Run pipeline as follows:

```
snakemake -n -d ../run_folder/ --cluster-status ../run_folder/hhu_hilbert.py --configfiles smk_config/smk_cfg_env-hhu.yml smk_config/kmervar_cfg.yml --profile environment/snakemake/cluster/hhu_pbs master
```

Configuration file `kmvar_cfg.yml` contains path to assemblies and k-mer size.

Current project folder:

```
/gpfs/project/projects/medbioinf/projects/kmvar
```

Input assemblies:

```
/gpfs/project/projects/medbioinf/projects/kmvar/assemblies
```

Adding more assemblies (FASTA format) to the above folder and (force) rerunning the pipeline will trigger an update for the k-mer counting.

The main output file of shared unique k-mers is located under:

```
run_folder/output/shared_kmers/uniq.k31.txt.gz
```
