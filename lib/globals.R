# Add any project specific configuration here.
add.config(
  apply.override = FALSE
)

# Add project specific configuration that can be overridden from load.project()
add.config(
  apply.override = TRUE,
  bioMart = 'data/mart_export.txt',
  diffRnaFile = 'reports/diffGenes.csv',
  lncRNA = c('lincRNA',
              'bidirectional_promoter_lncRNA',
              '3prime_overlapping_ncRNA',
              'macro_lncRNA',
              'antisense',
              'sense_overlapping',
              'sense_intronic'),
  PCGs = c('protein_coding',
            'IG_C_gene',
            'IG_J_gene',
            'TR_D_gene',
            'TR_V_gene',
            'TR_J_gene',
            'IG_D_gene',
            'IG_V_gene',
            'TR_C_gene')
)
