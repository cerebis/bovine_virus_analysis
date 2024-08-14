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

workflow clean_up {

  take:
  read_sets 

  main:
  Fastp(read_sets)
  remove_human(Fastp.out.reads)
  remove_host(remove_human.out)

  emit:
  remove_host.out
}

workflow {
    read_sets = Channel.fromPath(params.read_sets).splitCsv()
    clean_up(read_sets)
}
