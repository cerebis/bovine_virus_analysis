
params.asm_maxmem = '128 GB'
params.asm_cpus = 32
params.asm_queue = 'workq'


process MetaQuast {
    conda '/shared/homes/120274/miniconda3/envs/quast'
    publishDir params.outdir, mode: 'copy', saveAs: {fn -> "${name}/${asm_mode}spades/${fn}"}
    cpus 4
    memory '32 GB'
   
    input:
    tuple val(name), val(asm_mode), path(contigs)

    output:
    path('metaquast_out')

    """
    metaquast -t ${task.cpus} --max-ref-number 0 -o metaquast_out $contigs
    """
}

process Spades {
    scratch '/scratch/u120274/'
    conda '/shared/homes/120274/miniconda3/envs/spades'
    publishDir params.outdir, mode: 'copy', saveAs: {fn -> "${name}/${asm_mode}spades/"}
    cpus params.asm_cpus
    memory params.asm_maxmem
    queue params.asm_queue
    scratch '/scratch/u120274' 
    
    input:
    tuple val(name), path(r1), path(r2)
    val(asm_mode)

    output:
    tuple val(name), val(asm_mode), path('out/final.fasta'), emit: asm_seqs
    path('out'), emit: outdir

    script:
    if (asm_mode == "rna") {
        """
        spades.py --rna \
                  --threads ${task.cpus} \
                  -1 $r1 -2 $r2 \
                  -o out
        cd out
        mv transcripts.fasta final.fasta
        ln -s final.fasta transcripts.fasta
        """
    }
    else if (asm_mode == "meta") {
        """
        spades.py --meta \
                  --threads ${task.cpus} \
                  -1 $r1 -2 $r2 \
                  -o out
        cd out
        mv scaffolds.fasta final.fasta
        ln -s final.fasta scaffolds.fasta
        """
    }
}

process KmerSpades {
    conda '/shared/homes/120274/miniconda3/envs/spades'
    publishDir params.outdir, mode: 'copy', saveAs: {fn -> "${name}/${asm_mode}spades/"}
    cpus params.asm_cpus
    memory params.asm_maxmem
    queue params.asm_queue
    scratch '/scratch/u120274'

    input:
    tuple val(name), path(r1), path(r2)
    val(kmer_series)

    output:
    tuple val(name), val(asm_mode), path('out/final.fasta'), emit: asm_seqs
    path('out'), emit: outdir

    script:
    """
    spades.py --meta \
              -k $kmer_series \
              --threads ${task.cpus} \
              -1 $r1 -2 $r2 \
              -o out
    cd out
    mv scaffolds.fasta final.fasta
    ln -s final.fasta scaffolds.fasta
    """
}

process ExtendedSpades {

    conda '/shared/homes/120274/miniconda3/envs/spades'
    publishDir params.outdir, mode: 'copy', saveAs: {fn -> "${name}/${asm_mode}spades_extend/"}
    cpus 32
    memory '128 GB'

    input:
    tuple val(name), path(asm_work)
    val(asm_mode)
    val(kmers)

    output:
    tuple val(name), val(asm_mode), path('out/final.fasta'), emit: asm_seqs
    path('out'), emit: outdir
    

    script:
    if (asm_mode == "rna") {
        """
        spades.py -k $kmers \
                  --restart-from k73 \
                  -o $asm_work
        cd out
        mv transcripts.fasta final.fasta
        ln -s final.fasta transcripts.fasta
        """
    }
    else if (asm_mode == "meta") {
        """
        spades.py -k $kmers \
                  --restart-from k55 \
                  -o $asm_work
        cd out
        mv scaffolds.fasta final.fasta
        ln -s final.fasta scaffolds
        """
    }
}

workflow default_spades {

    take:
    read_sets
    asm_mode

    main:
    Spades(read_sets, asm_mode)

    emit:
    asm_seqs = Spades.out.asm_seqs
}

workflow kmer_spades {

    take:
    read_sets
    kmer_series

    main:
    KmerSpades(read_sets, kmer_series)

    emit:
    asm_seqs = KmerSpades.out.asm_seqs
}

workflow rna_spades {
    take:
    read_sets

    main:
    default_spades(read_sets, 'rna')

    emit:
    asm_seqs = default_spades.out.asm_seqs
}

workflow meta_qc {
    take:
    asm

    main:
    MetaQuast(asm)

    emit:
    MetaQuast.out
}

workflow meta_qc2 {
    take:
    asm

    main:
    MetaQuast(asm)

    emit:
    MetaQuast.out
}

workflow meta_spades {
    take:
    read_sets

    main:
    default_spades(read_sets, 'meta')

    emit:
    asm_seqs = default_spades.out.asm_seqs
}

workflow extended_spades {

    take:
    read_sets
    asm_mode
    extra_kmers

    main:
    Spades(read_sets)
    ExtendedSpades(Spades.out, asm_mode, extra_kmers)

    emit:
    asm_seqs = ExtendedSpades.out.asm_seqs
}
