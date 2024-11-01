process Sketch {
  conda '/shared/homes/120274/miniconda3/envs/sourmash'
  publishDir params.outdir, mode: 'copy', saveAs: {fn -> "$name/sourmash/$fn"}
  cpus 2
  memory '32 GB'

  input:
  tuple val(name), path(r1), path(r2)

  output:
  tuple val(name), path('sketch.sig.gz')
  
  """
  sourmash sketch dna -p k=31,abund $r1 $r2 -o sketch.sig.gz --name $name
  """
}

process Gather {
  conda '/shared/homes/120274/miniconda3/envs/sourmash'
  publishDir params.outdir, mode: 'copy', saveAs: {fn -> "$name/sourmash/$fn"}
  cpus 16
  memory '128 GB'

  input:
  tuple val(name), path(sketch)
  val(ref_db)
  val(min_len)

  output:
  tuple val(name), path('prefetch.zip'), path('gather.csv')

  """
  sourmash scripts fastgather \
        --cores ${task.cpus} \
        --threshold-bp $min_len \
        --output-prefetch prefetch.zip \
        --output-gather gather.csv \
        $sketch $ref_db
  """
}

process GatherExtra {
  conda '/shared/homes/120274/miniconda3/envs/sourmash'
  publishDir params.outdir, mode: 'copy', saveAs: {fn -> "$name/sourmash/$fn"}
  cpus 2
  memory '256 GB'

  input:
  tuple val(name), path(sketch), path(matches)
  val(ref_db)
  val(min_len)

  output:
  tuple val(name), path('gather_extra.csv')
  path('matches_extra.zip') 

  """
  sourmash gather --threshold-bp $min_len $sketch $ref_db --save-matches matches_extra.zip
  sourmash gather --threshold-bp $min_len $sketch $matches matches_extra.zip -o gather_extra.csv
  """
}

process TaxProfile {
  conda '/shared/homes/120274/miniconda3/envs/sourmash'
  publishDir params.outdir, mode: 'copy', saveAs: {fn -> "$name/sourmash/$fn"}
  cpus 2
  memory '64 GB'

  input:
  tuple val(name), path(matches), path(gather)
  val(tax_db)
  val(rank)

  output:
  path('taxprof.*')

  """
  sourmash tax metagenome \
        -F human csv_summary lineage_summary krona kreport \
        -g $gather \
        -t $tax_db \
        -r $rank \
        -o taxprof
  """
}

process TaxProfileExtra {
  conda '/shared/homes/120274/miniconda3/envs/sourmash'
  publishDir params.outdir, mode: 'copy', saveAs: {fn -> "$name/sourmash/$fn"}
  cpus 2
  memory '64 GB'

  input:
  tuple val(name), path(gather)
  val(tax_db)
  val(rank)

  output:
  path('taxprof_extra.*')

  """
  sourmash tax metagenome \
        -F human csv_summary lineage_summary krona kreport \
        -g $gather \
        -t $tax_db \
        -r $rank \
        -o taxprof_extra
  """
}

process CAT {
  conda '/shared/homes/120274/miniconda3/envs/catbat'
  scratch '/scratch/u120274'
  publishDir params.outdir, mode: 'copy', saveAs: {fn -> "$name/cat/$asm_run/$fn"}
  cpus 32
  memory '250 GB'

  input:
  tuple val(name), val(asm_run), path(contigs)
  val(cat_db)
  val(cat_tax)

  output:
  tuple val(name), val(asm_run), path('CAT*')

  """
  /shared/homes/120274/git/CAT_pack/CAT_pack/CAT_pack contigs \
    -n ${task.cpus} --force --compress \
    -d $cat_db \
    -t $cat_tax \
    -c $contigs \
    -o CAT
  """
}

process CAT_addnames {
  conda '/shared/homes/120274/miniconda3/envs/catbat'
  publishDir params.outdir, mode: 'copy', saveAs: {fn -> "$name/cat/$asm_run/$fn"}
  cpus 2
  memory '32 GB'

  input:
  tuple val(name), val(asm_run), path(classification)
  val(cat_tax)

  output:
  tuple val(name), val(asm_run), path('*.names.tsv')

  shell:
  '''
  for fn in CAT.contig2classification.txt CAT.ORF2LCA.txt;
  do
      /shared/homes/120274/git/CAT_pack/CAT_pack/CAT_pack add_names \
        --force --exclude_scores \
        -t !{cat_tax} \
        -i $fn \
        -o "$(basename $fn .txt).names.tsv"
  done
  '''
}

workflow contig_assign {
    take:
    contig_sets
    cat_db
    cat_tax

    main:
    CAT(contig_sets, cat_db, cat_tax)
    CAT_addnames(CAT.out, cat_tax)

    emit:
    CAT.out
    CAT_addnames.out
}

workflow reads_profile {
  take:
  read_sets
  ref_db
  tax_db
  min_len
  rank

  main:
  Sketch(read_sets)
  Gather(Sketch.out, ref_db, min_len)
  TaxProfile(Gather.out, tax_db, rank)

  emit:
  TaxProfile.out
}



workflow {

    sketches = Channel.fromPath(params.sketches).splitCsv()

    GatherExtra(sketches, 
                Channel.value(params.reference_db))

    TaxProfileExtra(GatherExtra.out[0], 
                    Channel.value(params.taxonomy_db), 
                    Channel.value(params.rank))


//  contig_sets = Channel.fromPath(params.contig_sets).splitCsv()
//  contig_assign(contig_sets, Channel.fromPath(params.cat_db), Channel.fromPath(params.cat_tax))

//  cat_sets = Channel.fromPath(params.cat_sets).splitCsv()
//  CAT_addnames(cat_sets, Channel.fromPath(params.cat_tax))
}
