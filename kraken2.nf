process Kraken2 {
    publishDir params.outdir, mode: 'copy', saveAs: {fn -> "${run}/${fn}"}
    cpus 10
    conda '/shared/homes/120274/miniconda3/envs/kraken2'
    executor 'local'
    scratch '/scratch/u120274'
   
    input:
    tuple val(run), path(r1), path(r2)
    val(db)

    output:
    tuple path('output.krk'), path('report.tsv')

    """
    kraken2 --memory-mapping \
            --threads ${task.cpus} \
            --paired \
            --output output.krk \
            --report report.tsv \
            --db $db \
            $r1 $r2
    """
}

workflow {
    read_sets = Channel.fromPath(params.runs).splitCsv()
    Kraken2(read_sets, Channel.value(params.db))
}
