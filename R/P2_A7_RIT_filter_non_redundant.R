### P2_A7 --------------------------------------------------------------------

## This script does:
## Removing redundant INSD records for which there is already a RefSeq record
## Fixing taxonomy labels where missing.

### NOTES ----------------------------------------------------------------

## Nucleotide summaries are missing too many assembly GIs and accessions to use those for filtering.

## Checking:
# I checked all the records by hand and made a giant list of duplicated records where both INSD and RefSeq had records for the same genome for a given strain.

## Tried identify any assemblies represented in insd and refseq where refseq is missing a rit.
# INSD sequences for taxa not represented in Refseq subset.
# tax_id 169430
# 997885, RIT_14, Bacteroides ovatus CL02T12C04

## Other notes:
# A back to back dimer RIT
# https://www.ncbi.nlm.nih.gov/nuccore/UGYX01000003.1?report=genbank&log$=seqview&from=734712&to=742000

# BROKEN ASSEMBLY: ---- strain names for assemblies
# Rit removed from Refseq assembly - present in INSD assembly
# Could just be missing corresponding RefSeq records though - not necessarily missing from the RefSeq assembly:
# "Bacteroides finegoldii DSM 17565" ids:
# 238852439 240114316 242353479
# 'Bacteroides sp. 4_1_36' refseq is missing insd contig:
# 316902394 318103433 319430555
# Caballeronia sordidicola strain ES_PA-B12 Ga0193686_11 - not in refseq?
# 1534957005 1534954659 1537750292 1534944218
# Johnsonella ignava ATCC 51276,
# Marinobacter sp. ELB17 1101232001194,
# Prevotella buccae ATCC 33574 contig00043 is missing from RefSeq?
# Acidovorax facilis isolate 7, whole genome shotgun sequence - has 2 different assemblies, with different RITs.
# Cupriavidus taiwanensis isolate Cupriavidus taiwanensis mpp 1.1 plasmid CBM2626_p
# has a plasmid RIT that is absent in assembly NZ_LT976998.1; 

# RECENT MOBILITY NOTES -----
# multiple copies of identical element. 
# Croceicoccus marinus strain OT19 chromosome and plasmid plas1 (3 copies + another RIT)
# Rhizobium leguminosarum strain ATCC 14479 chromosome and plasmids x 6.
# Sphingobium herbicidovorans strain MH - chromosome and plasmids
# Sphingobium sp. TKS 

# Setup =================================================================

library(tidyverse)
library(glue)
library(lubridate)

for_viewing <- function(rits){
  rits |> 
    select(-contains('dna')) |> 
    relocate(tax_id, rit_id, organism, strain, sourcedb, title,
             slen, rit_start, rit_stop) |> 
    View()
}


# Id lists =================================================================

# This is a list of redundant ids to remove, already have a different assembly.
# Prefer refseq to insd where similar record is given by both; prefer scaffold to contig-level assembly
# I noticed that sometimes where insd and refseq have same assemblies - insd has extra contigs that also contain rits - have left these in place.
# This could just be because I didn't get ALL the possible genbank records

redundant_records <- c(
  '400374209', '654778780', '918160129', '355544533', '1044497026',
  '1549570730', '1549652570', '1549634786', '1549695024', '786805399',
  '1622219689', '685434132', '1631997304', '1778136453', '1097499250',
  '1593380907', '1254831335', '645571062', '1440787348', '742252786',
  '524049784', '1778136473', '1778136477', '1198728346', '1240584203',
  '1240583949', '1240584203', '1240583372', '1240583949', '1573062409',
  '1436261838', '700235774', '1369277669', '1047883317', '1098213502',
  '1109830987', '765543252', '1132730147', '1297911868', '1703336581',
  '1780515536', '512186279', '805451489', '805426424', '805451498',
  '1573062409', '735022113', '1418835274', '698841544', '1596424614',
  '1716017357', '1716010656', '1515626905', '1515608535', '1515607855', 
  '1485368893','1485366470', '1496378801', '1485377885', '1485366962',
  '1596433374', '1216427170', '1216429375', '724421287', '751285260', 
  '1618902489', '1442853557', '1443956302', '1442886711', '1442853162',
  '1097189591', '1377955386', '1393507009',  '962231498', '1540051836', 
  '1405694183', '563003182', '1622715784', '1173527474', '1208465567',
  '735022113', '507806472', '512222473', '805426424', '1373130578',
  '1197485940', '1197459102', '1736577639', '1736577573', '1736591669', 
  '1488751460', '1596446976', '1485352105',  '1485361253', '1596216906',
  '1596439722', '1596208487', '1488751460', '928763215', '928737052', 
  '1596318712', '397151626', '392661434', '913264544', '805435520', 
  '805426915', '397151634', '1596096393', '1596098821', '1596104534', 
  '939504136', '752200082', '998167973', '1269187544', '1436254977',
  '595611723', '945493622', '1087612761', '1087613108', '1597161555', 
  '1446762304', '1446761670', '1485450145', '1354608460', '1024316281',
  '507807538', '13488050', '1621679256', '646767384',
  '1252079409', '1252085149', '1252078968', '1356327944', '392663357',
  '1442852763', '1443253504', '1198804558', '1500099130', '1198804472',
  '1198444759','1716021005','1485372803','1515632920','1780515536',
  '1442852373','1418835274','698178797','389931344','1596208487', 
  '805430678','1280122714','BL6665CI2N2','512222474', '978362936',
  '1252082629', '169430', '93218', '1081889900', '1032499711', '1088980949',
  '363417490', '946165513', '940794514', '1703348145', '1703348283', 
  '1450641438', '1046102153', '1215060541', '221738762', '1198774940',
  '1382859007', '1383481127', '1004898246', '1470412422', '1470401030',
  '1745926388', '1469836909' ,'1469844101', '1470686740', '1467104658',
  '397151430', '397151492', '596159047', '596170443', '598872967', '598897908',
  '596068870', '596170915', '595938801', '596273949', '596107824', '595946098',
  '595952780', '1290342027', '1720182277', '1764098768', '1720189434', 
  '1471907645', '1047137029', '1047135602', '1720208333', '1046946578',
  '1046966569', '1046943050', '1046966570', '1087270618', '318103433',
  '316902394', '1484146896', '1484315991' ,'1198643640', '1198654165',
  '1198492777', '1110476557', '1470361372', '382984480', '1470138912', 
  '1469878467', '1762752354', '1762762889', '1762762945' ,'1762751577',
  '649533298', '1470967598', '1198794869', '1198793338', '1765288693', '1765260936',
  '1765269603' ,'1765263374', '1765295756', '1765240997', '1765233681', '1765248576',
  '1765253103', '1765227581', '1765233347', '1765247803', '1765224780', '1765215305',
  '1570814390', '1099534514', '1466548916', '1198502171', '1569619122', '983695177',
  '983693665', '1562167019', '1759753624', '1533368274', '1699440712', '1177748598',
  '1525273214', '1525121968', '636800476', '1252060762', '831716647', '1534954659' ,
  '1534944218', '1588246170', '1822154115', '1815527696', '1815531171', 
  '752651577', '752659692', '524724042', '1110268723', '1131990612', '972591476',
  '1108333464', '1116070613', '671228823', '1146522581', '1159623998', '1098328891',
  '1593032722', '1418276311', '530294384', '530293375', '530294384', '530293375', 
  '1148971634', '1639957942' ,'1639957071', '1639956401', '1639955835', '1446356571',
  '1241768718', '1153801033', '1063644618', '945697828', '945731194', '946324979',
  '1360186106', '1198689136', '1472090583', '512600894', '1198755451' ,'1198645385',
  '1547682791', '1387694158', '1691776874', '1407648773', '1407648769', '1596940895', 
  '1596940888', '1086611919', '1289242642', '1465780481', '357379774', '1198538775',
  '363537050', '962209145', '962247357', '1437693291' ,'999165770', '962276359', 
  '962274924', '1437550231', '1098092310', '126626772', '126626409', '1481335117',
  '1692135553', '1692147915' ,'1692130005', '1692162480' , '1692151059' , '1692165367',
  '1692140722', '1692113932' , '1692076492', '1692079481' , '1692100917', '1692057266',
  '1692009440' , '1547976866', '1548090616', '1614958144', '1614949718', '1743057145',
  '1743055600', '1743544331', '1743544320', '1099823037', '1028357265', '929665776',
  '1175032771', '1399513828', '808402906', '1470672900', '805451493', '1226586239',
  '1466733420', '1238362427', '1722658860', '1431882981', '1431875858', '1440788695',
  '1597162718', '1588161165', '1397780655', '1405754217', '1488931450' ,'1574044243',
  '1110540023', '1110537708', '1120902308', '1120902418', '1405831238', '1397906827',
  '529584578', '1680504601', '1087359298', '1715802498', '1086648911', '945497149', 
  '945594613', '946022763', '1212230502', '1471732630', '1471713537', '1470627261',
  '1470620688', '1469566809', '1471910516', '1471583085', '1470970509', '1708670197',
  '1708675796', '1466603683', '317575648', '315249706', '288336044', '1437513354',
  '1545075283', '1097302153', '1823128862', '1586509317', '1586479078', '835569958',
  '835569972', '1248503826', '1487678657', '1254765167' ,'1669254997', '1099798388', 
  '1485389667', '1580603362', '1470686166', '1731234072', '1699357987', '1472023606',
  '1472789073' , '1472765508', '1109644811', '1578214595', '1589256909', '1488048787',
  '1564781530', '1011990542', '1011990577', '1011990580', '674786029',
  '1779023398', '1779023566', '1269808047', '1125936988', '1537305031', 
  '1537247517', '1537218261', '1537237467', '1537661137' , '1537644841', '1537222108',
  '1537658515', '1537671529', '1537647300', '1537630491' ,'1537166756', '1537632099',
  '1537643148', '1537657706' , '1537619749', '1537618519', '1537610293', '1537599012',
  '1537602337', '1537249002', '1537271918', '1537268662', '1537192810' , '1537175087', 
  '1537182942', '1537289869', '1537281862' , '1537295193', '1537441230', '1537251846',
  '1537229864', '1537167927','1537456472', '1537456472','1537452925', '1537453651', 
  '1537462788', '1537476592', '1537467426', '1537487018', '1537490750' , '1537301764',
  '1537509386', '1537525587', '1537483753', '1537514509', '1537521166', '1537510319',
  '1537538998', '1537534600', '1537538420', '1537529883', '1537163295', '1537554636',
  '1537558745', '1537567479', '1537547964', '1537573649', '1537558806', '1537570678',
  '1537587365', '1537585139', '1537593289', '1537180379', '1537580224', '1537607698',
  '1537602337', '1537599012', '1537610293', '1537618519' , '1537619749', '1537657706',
  '1537643148', '1537632099', '1537166756', '1537630491', '1537647300', '1537671529',
  '1537658515', '1537222108', '1537644841', '1537661137', '1537237467', '1537218261',
  '1537247517','1537305031', '1730169111', '1631418225', '1253817562', '1588658013', 
  '1072721238', '1072722420' , '821158070', '1691230684', '1024470410', '1024466825',
  '1484222692', '995259649','962388704', '962340824', '962321783' , '962378348', 
  '962337292', '962354456', '962322963', '962347846', '962316371', '962356313',
  '1417744410', '1133385480', '1537203256', '946135381', '945815692' , '469481440', 
  '1700681327', '1537469923', '384039839', '1174162640', '1019513169', '1019513283',
  '1019514939', '1019395847', '1253857670', '25168258', 
  '356692854', '942698728', '1777945395', '392627216', '966499388', '966499399', 
  '1213754548', '821179869', '754936391', '918697158', '1280123531', '507781804'
)

# Was starting to add strain info where present in title but not in proper column for strain from the taxonomy records. Gave up. Too many to fix by hand.
strains <- tribble(
  ~nuc_id, ~strain,
  '1899015', 'MES1',
  '1262767', 'CAG:290',
  '1371401794', 'AMDSBA3'
)



## 1 Data ---------------------------------------------

# 2,262 records from rit_finder
rit_elements <- read_rds('./results/rit_elements.rds') |> 
  select(-trio_id)
glimpse(rit_elements)

## Keep only the refseq records for duplicated (seq, taxon) pairs.
# 1,558  RITs after filtering redundant records
rit_elements <- rit_elements |> 
  filter(!nuc_id %in% c(redundant_records)) 
glimpse(rit_elements)

# distinct_rits -> obs grouped by RIT seq. from start p1 to end p3 and taxon id
# rits grouped by exact dna sequence over coding region
distinct_rits <- read_rds('results/distinct_RITs.rds')
glimpse(distinct_rits)

# records summary data to join to filtered rits
nuc_summary <- read_rds('data/CDD/nuc_summary_fixed.rds') |> 
  select(-c(replacedby, segsetsize, tech, moltype, oslt, status, strand,
            subtype, subname, statistics, properties, accessionversion,
            caption, geneticcode, extra, idgiclass, comment)) |> 
  mutate(slen = as.numeric(slen),
         across(contains('date'), ~lubridate::ymd(.x))) 
glimpse(nuc_summary)


# 666 distinct taxa
rit_elements |> dplyr::count(tax_id) |> nrow()
# 1086 nr records
rit_elements |> dplyr::count(nuc_id) |> nrow()
# 916 unique RIT sequences from p1-p3
rit_elements |> dplyr::count(rit_dna) |> nrow()

## 2 Filter -------------------------------------------

# unnest rit elements' occurrences in all records
# keep non-redundant records
# join nucleotide summary data
# update rit_id using the unique ids from rit_unique_seq_id
rits_with_nuc_details <-  distinct_rits |> 
  unnest(cols = c(rit_occurrences)) |> 
  left_join(nuc_summary, by = 'nuc_id') |> 
  mutate(seq_type = str_extract_all(tolower(nuc_name),
                                    'contig|scaffold|complete genome'),
         strain = coalesce(strain.y, strain.x))  |> 
  mutate(rit_id = rit_unique_seq_id) |> 
  select(-c(rit_unique_seq_id, trio_id, tax_genetic_code, 
            superkingdom, strain.x, strain.y))
glimpse(rits_with_nuc_details)

# these are the 704 records that were removed...
redundant_rec <-  rits_with_nuc_details |> 
  filter(nuc_id %in% c(redundant_records))
redundant_rec
write_rds(redundant_rec, 'data/redundant_records_rits.rds')

# 17 that are redundant assemblies...
# redundant_asm <-  rits_with_nuc_details |> 
#   filter(nuc_id %in% c(redundant_assemblies))
# redundant_asm
# write_rds(redundant_asm, 'data/redundant_assembly_rits.rds')

# 1558 rit results across 1086 nr records
nr_rits <- rits_with_nuc_details |> 
  filter(!nuc_id %in% c(redundant_records)) 
nr_rits |> dplyr::count(nuc_id) |> nrow()

## visually inspect sequence records included in data w list of RITs
## for_viewing(nr_rits)


# fix missing taxonomy labels
rits <- nr_rits |> 
  mutate(
    phylum = case_when(
      str_detect(phylum, 'Spirochaetes') ~ 'Spirochaetes',
      str_detect(phylum, 'Acidobact') ~ 'Acidobacteria',
      str_detect(phylum, 'Bacteria candidate phyla') ~ 'Candidate',
      is.na(phylum) ~ 'unclassified',
      TRUE ~ phylum
      ),
    phylum = str_remove_all(phylum, ' group| subdivisions'),
    phylum = str_to_title(phylum),
    class = ifelse(is.na(class), glue('unclassified {phylum}'), class),
    class = str_replace(class, 'Betaproteobacteria incertae sedis',
                        'unclassified Betaproteobacteria')
    ) 

rits |> 
  select(tax_id, clade, phylum, class, order, family, genus, species) |> 
  distinct() |> 
  dplyr::count(clade, phylum, class, order, family, genus, species) |> 
  arrange(desc(n)) |> 
  print(n=100)

glimpse(rits)

# write to file for analysis: finding active RITs, finding IRs
write_rds(rits, 'results/non_redundant_rits.rds')

beepr::beep()

# View(
#   nr_rits |>
#     group_by(tax_id, tax_name, strain, nuc_id, title, sourcedb,
#              slen, completeness, seq_type, assemblygi, assemblyacc, projectid) |>
#     summarize(rits = paste0(rit_id, collapse = ', ')) |>
#     distinct() |>
#     arrange(tax_id, sourcedb)
# )
# 
# 
