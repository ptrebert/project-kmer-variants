import os
import sys

localrules: master

INPUT_FASTA_PATH = config['input_fasta_path']
KMER_SIZE = config['kmer_size']

ASSEMBLIES = []
for item in os.listdir(INPUT_FASTA_PATH):
    if not item.endswith('.fasta'):
        sys.stderr.write('\nSkipping over {}...'.format(item))
    ASSEMBLIES.append(item.rsplit('.', 1)[0])


rule count_kmers_in_assembly:
    """
    Create k-mer db
    and also a dump highly frequent k-mers
    (could be needed for winnowmap alignment to that assembly)
    """
    input:
        fasta = os.path.join(INPUT_FASTA_PATH, '{assembly}.fasta')
    output:
        kmer_db = directory('output/kmer_db/{assembly}.k{kmer}.db/'),
        rep_kmer = 'output/kmer_freq/{assembly}.k{kmer}.rep-grt09998.txt'
    benchmark:
        'rsrc/output/kmer_db/{assembly}.k{kmer}.count.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_medium']
    resources:
        mem_total_mb = lambda wildcards, attempt: 32768 * attempt,
        mem_total_gb = lambda wildcards, attempt: 32 * attempt,
        runtime_hrs = lambda wildcards, attempt: attempt
    shell:
        'meryl count k={wildcards.kmer} threads={threads} memory={resources.mem_total_gb} output {output.kmer_db} {input.fasta} '
        '&& '
        'meryl print greater-than distinct=0.9998 {output.kmer_db} > {output.rep_kmer}'


def assemble_meryl_intersect(kmer_dbs):

    if isinstance(kmer_dbs, str):
        # single DB would fail anyway
        raise ValueError('Single k-mer DB detected, meryl intersect needs at least two inputs: {}'.format(kmer_dbs))
    template = '[equal-to 1 {}] '
    call = ''
    for db in kmer_dbs:
        call += template.format(db)
    return call


rule determine_shared_kmers:
    input:
        expand(
            'output/kmer_db/{assembly}.k{{kmer}}.db/',
            assembly=ASSEMBLIES,
        )
    output:
        db = directory('output/shared_kmers/uniq.k{kmer}.db/'),
        listing = 'output/shared_kmers/uniq.k{kmer}.txt'
    log:
        'log/output/shared_kmers/uniq.k{kmer}.meryl.log'
    benchmark:
        'rsrc/output/shared_kmers/uniq.k{kmer}.meryl.rsrc'
    conda:
        '../../environment/conda/conda_biotools.yml'
    threads: config['num_cpu_medium']
    resources:
        mem_total_mb = lambda wildcards, attempt: 4096 * attempt,
        mem_total_gb = lambda wildcards, attempt: 4 * attempt,
        runtime_hrs = lambda wildcards, attempt: 8 * attempt
    params:
        meryl_call = lambda wildcards, input: assemble_meryl_intersect(input)
    shell:
        'meryl threads={threads} memory={resources.mem_total_gb} '
        'intersect '
        '{params.meryl_call} '
        'output {output.db} &> {log} '
        '&& '
        'meryl print {output.db} > {output.listing}'


rule master:
    input:
        expand(rules.determine_shared_kmers.output.listing, kmer=KMER_SIZE)
