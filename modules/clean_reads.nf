process FastQC {
    conda '/shared/homes/120274/miniconda3/envs/fastqc'
    publishDir params.outdir, mode: 'copy', saveAs: {fn -> "$name/qc/reads/$stage/$fn"}
    cpus 8
    memory '16 GB'
    
    input:
    tuple val(name), path(r1), path(r2)
    val(stage)

    output:
    path('fastqc_out')

    """
    mkdir fastqc_out
    fastqc -t ${task.cpus} -o fastqc_out $r1 $r2
    """
}

process Fastp {
  conda '/shared/homes/120274/miniconda3/envs/fastp'
  publishDir params.outdir, mode: 'copy', saveAs: {fn -> "$name/fastp/$fn"}
  cpus 16
  memory '32 GB'

  input:
  tuple val(name), path(r1), path(r2)

  output:
  tuple val(name), path('out/*1.fastq.gz'), path('out/*2.fastq.gz'), emit: reads
  path('out/report.html')

  """
  mkdir -p out
  fastp -w ${task.cpus} --detect_adapter_for_pe --dedup \
        -f 1 -F 1 \
        -i $r1 -I $r2 \
        -o out/${r1.fileName} -O out/${r2.fileName} \
        -h out/report.html
  """
}

process Hostile {
  conda '/shared/homes/120274/miniconda3/envs/hostile'
  publishDir params.outdir, mode: 'copy', saveAs: {fn -> "$name/hostile/$refname/$fn"}
  cpus 16
  memory '64 GB'

  input:
  tuple val(name), path(r1), path(r2)
  val(refname)
  val(index)

  output:
  tuple val(name), path('out/*1.fastq.gz'), path('out/*2.fastq.gz'), emit: reads
  path('hostile.log')

  shell:
  '''
  hostile clean --fastq1 !{r1} --fastq2 !{r2} \
                --index !{index} \
                --threads !{task.cpus-3} \
                --out-dir out > hostile.log
  '''
}

workflow remove_human {

  take:
  read_sets

  main:
  Hostile(read_sets, 'stage1', Channel.value(params.human_index))

  emit:
  Hostile.out.reads
}

workflow remove_host {

  take:
  read_sets

  main:
  Hostile(read_sets, 'stage2', Channel.value(params.host_index))

  emit:
  Hostile.out.reads
}

workflow trim {
    take:
    read_sets

    main:
    Fastp(read_sets)

    emit:
    Fastp.out.reads
}

workflow fastqc_raw {
    take:
    read_sets

    main:
    FastQC(read_sets, Channel.value('raw'))
}

workflow fastqc_clean {
    take:
    read_sets

    main:
    FastQC(read_sets, Channel.value('clean'))
}

workflow fastqc_nohuman {
    take:
    read_sets

    main:
    FastQC(read_sets, Channel.value('no_human'))
}

workflow fastqc_nocow {
    take:
    read_sets

    main:
    FastQC(read_sets, Channel.value('no_cow'))
}

workflow clean_up {

  take:
  read_sets 

  main:

      fastqc_raw(read_sets)

  Fastp(read_sets)

      fastqc_clean(Fastp.out.reads)

  remove_human(Fastp.out.reads)

      fastqc_nohuman(remove_human.out)

  remove_host(remove_human.out)

     fastqc_nocow(remove_host.out)

  emit:
  remove_host.out
}

workflow {
    read_sets = Channel.fromPath(params.read_sets).splitCsv()
    clean_up(read_sets)
}
