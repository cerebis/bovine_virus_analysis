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
  cpus 2
  memory '256 GB'

  input:
  tuple val(name), path(sketch)

  output:
  tuple val(name), path('matches.zip'), path('gather.csv')

  """
  sourmash gather --threshold-bp $params.minlen $sketch ${params.reference_db} --save-matches matches.zip
  sourmash gather --threshold-bp $params.minlen $sketch matches.zip -o gather.csv
  """
}

process GatherExtra {
  conda '/shared/homes/120274/miniconda3/envs/sourmash'
  publishDir params.outdir, mode: 'copy', saveAs: {fn -> "$name/sourmash/$fn"}
  cpus 2
  memory '256 GB'

  input:
  tuple val(name), path(sketch), path(matches)
  val(refdb)

  output:
  tuple val(name), path('gather_extra.csv')
  path('matches_extra.zip') 

  """
  sourmash gather --threshold-bp $params.minlen $sketch $refdb --save-matches matches_extra.zip
  sourmash gather --threshold-bp $params.minlen $sketch $matches matches_extra.zip -o gather_extra.csv
  """
}

process TaxProfile {
  conda '/shared/homes/120274/miniconda3/envs/sourmash'
  publishDir params.outdir, mode: 'copy', saveAs: {fn -> "$name/sourmash/$fn"}
  cpus 2
  memory '64 GB'

  input:
  tuple val(name), path(matches), path(gather)

  output:
  path('taxprof.*')

  """
  sourmash tax metagenome -F human csv_summary lineage_summary krona kreport \
     -g $gather -t ${params.taxonomy_db} -r ${params.rank} -o taxprof
  """
}

process TaxProfileExtra {
  conda '/shared/homes/120274/miniconda3/envs/sourmash'
  publishDir params.outdir, mode: 'copy', saveAs: {fn -> "$name/sourmash/$fn"}
  cpus 2
  memory '64 GB'

  input:
  tuple val(name), path(gather)
  val(taxdb)
  val(rank)

  output:
  path('taxprof_extra.*')

  """
  sourmash tax metagenome \
        -F human csv_summary lineage_summary krona kreport \
        -g $gather \
        -t $taxdb \
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
  each path(cat_db)
  each path(cat_tax)

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
  each path(cat_tax)

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

  main:
  Sketch(read_sets)
  Gather(Sketch.out)
  TaxProfile(Gather.out)

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
