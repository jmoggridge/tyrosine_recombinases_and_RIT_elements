queries <- list(
  'Serine recombinase' = "serine recombinase",
  'Holliday junction resolvase' = 'Holliday junction resolvase',
  'Topoisomerase' = 'topoisomerase',
  'DDE transposase' = "DDE transposase",
  'DNA helicase recQ' = 'DNA helicase recQ',
  'Thermonuclease' = 'thermonuclease',
  'Ribonuclease' = 'ribonuclease',
  'Helicase' = 'helicase',
  'Histone' = 'histone',
  'Argonaute' = 'argonaute',
  'Restriction enz' = 'restriction endonuclease',
  'MbeA' = 'DNA relaxase mbeA domain protein',
  'RpnA' = 'Recombination-promoting nuclease RpnA',
  'RpnB' = 'Recombination-promoting nuclease RpnB',
  'RpnC' = 'recombination-promoting nuclease RpnC',
  'GamL' = 'host nuclease inhibitor GamL',
  'Nuclease SbcCD subunit D' = 'Nuclease SbcCD subunit D',
  'NikB' = 'NikB',
  'Exodeoxyribonuclease' = 'Exodeoxyribonuclease',
  'RecB' = 'exodeoxyribonuclease V subunit beta',
  'RecC' = 'exodeoxyribonuclease V subunit gamma',
  'RecD' = 'exodeoxyribonuclease V subunit alpha',
  'RecJ' = 'single-stranded-DNA-specific exonuclease RecJ',
  'RecN' = 'DNA repair protein RecN',
  'RecR' = 'recombination protein RecR',
  'Exonuclease V' = 'Exonuclease V',
  'IS607' = 'IS607 family transposase',
  'DnaX' = 'DNA polymerase III subunit gamma/tau',
  'Gntr' = 'gntr family transcriptional regulator',
  'Adenylosuccinate lyase' = 'Adenylosuccinate lyase',
  'IhfA' = 'integration host factor subunit alpha',
  'IhfB' = 'Integration host factor subunit beta',
  'RadA' = 'DNA repair protein RadA'
) |> 
  enframe(value = 'Search expression')

queries |> unnest(`Search expression`) |>  select(`Search expression`) |> 
  gt::gt() |> 
  gt::tab_options(table.font.size = 11, data_row.padding = 1)

queries |> unnest(`Search expression`) |>  pull(`Search expression`) |> 
  paste0(collapse = '; ')


