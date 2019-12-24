# R script to download selected samples
# Copy code and run on a local machine to initiate download

## set working directory
setwd('/home/locutus-of-borg/Documents/brent_lab/archs4_human')

# Check for dependencies and install if missing
packages <- c("rhdf5", "dplyr")
## include custom dependences functions here
#BiocManager::install("GEOquery")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  print("Install required packages")
  source("https://bioconductor.org/biocLite.R")
  biocLite("rhdf5")
}
library("rhdf5")
library("tools")
library("dplyr")
library("GEOquery")

destination_file = "human_matrix_download.h5"
extracted_expression_file = "human_possible_txfactor_datasets.csv"
url = "https://s3.amazonaws.com/mssm-seq-matrix/human_matrix.h5"

# Check if gene expression file was already downloaded and check integrity, if not in current directory download file form repository
if(!file.exists(destination_file)){
  print("Downloading compressed gene expression matrix.")
  download.file(url, destination_file, quiet = FALSE)
} else{
  print("Verifying file integrity...")
  checksum = md5sum(destination_file)
  
  if(destination_file == "human_matrix_download.h5"){
    # human checksum (checksum is for latest version of ARCHS4 data)
    correct_checksum = "34197866d7841cc4fb31e09195faa150"
  } else{
    # mouse checksum (checksum is for latest version of ARCHS4 data)
    correct_checksum = "55441d1af9da82c6f3d368c8fa554d42"
  }
  
  if(checksum != correct_checksum){
    print("Existing file looks corrupted or is out of date. Downloading compressed gene expression matrix again.")
    download.file(url, destination_file, quiet = FALSE)
  } else{
    print("Latest ARCHS4 file already exists.")
  }
}

checksum = md5sum(destination_file)
if(destination_file == "human_matrix_download.h5"){
  # human checksum (checksum is for latest version of ARCHS4 data)
  correct_checksum = "34197866d7841cc4fb31e09195faa150"
} else{
  # mouse checksum (checksum is for latest version of ARCHS4 data)
  correct_checksum = "55441d1af9da82c6f3d368c8fa554d42"
}

if(checksum != correct_checksum){
  print("File download ran into problems. Please try to download again. The files are also available for manual download at http://amp.pharm.mssm.edu/archs4/download.html.")
} else{
  # Selected samples to be extracted
  cas9 = c("GSM1145141","GSM1145142","GSM1145139","GSM1145138","GSM1145140","GSM2859897","GSM2136847","GSM2136848","GSM2136849","GSM2136850","GSM2136851","GSM2136852","GSM2136853","GSM2136854","GSM2136855","GSM2136856","GSM2136857","GSM2136858","GSM2136859","GSM2136860","GSM2136861","GSM2136862","GSM2136863","GSM2136864","GSM2136865","GSM2136866","GSM2136867","GSM2136868","GSM2136869","GSM2136870","GSM2136871",
           "GSM2136872","GSM2136873","GSM2136874","GSM2136875","GSM2136876","GSM2136877","GSM2136878","GSM2136879","GSM2248003","GSM2248004","GSM2248005","GSM2248006","GSM2248007","GSM2248008","GSM2248009","GSM2248010","GSM2248011","GSM2248012","GSM2248013","GSM2248014","GSM2248015","GSM2248016","GSM2248017","GSM2248018","GSM2248019","GSM2248020","GSM2248021","GSM2248022","GSM2248023","GSM3004226",
           "GSM3004227","GSM3004228","GSM3004229","GSM2572604","GSM2572605","GSM2572606","GSM2572607","GSM2572608","GSM2572609")
  crispr = c("GSM2520578","GSM3004226","GSM3004227","GSM3004228","GSM3004229")
  
  induction = c("GSM1829820","GSM1482170","GSM1829828","GSM1829830","GSM1482172","GSM1482171","GSM1482169","GSM1829826","GSM1829832","GSM1482166","GSM1482167","GSM1829822","GSM1829818","GSM1829809","GSM1482173","GSM1482174","GSM1829814","GSM1829824","GSM1829812","GSM1482168","GSM1829816","GSM2563696","GSM2563697","GSM2359781","GSM2359782","GSM2359783","GSM2359784","GSM2359785","GSM2359786","GSM2359787","GSM2359788",
                "GSM2359789")
  knockdown = c("GSM741170","GSM1254458","GSM1173504","GSM1254459","GSM1254457","GSM1294148","GSM2033072","GSM1294150","GSM1668429","GSM1508025","GSM1838574","GSM1838576","GSM1508029","GSM1668428","GSM1294151","GSM1668426","GSM1275184","GSM1838578","GSM1668432","GSM1668431","GSM1548041","GSM1548045","GSM1323580","GSM1267861","GSM1508003","GSM1668427","GSM1548040","GSM2056803","GSM1316405","GSM2056797","GSM1508031",
                "GSM1548038","GSM1294146","GSM2056802","GSM1548046","GSM1294149","GSM1548037","GSM2056800","GSM2056801","GSM1508001","GSM1668430","GSM1267860","GSM1838577","GSM1507999","GSM2056798","GSM1508026","GSM1838575","GSM2056799","GSM1684872","GSM1316404","GSM1508042","GSM1508018","GSM1508004","GSM1548044","GSM1267849","GSM1548042","GSM1508006","GSM1316406","GSM2056795","GSM1668433","GSM1294152",
                "GSM1508012","GSM1323581","GSM1508041","GSM1508027","GSM2056794","GSM1267855","GSM2056796","GSM1294147","GSM1548039","GSM1508040","GSM1294145","GSM1267848","GSM1508014","GSM2029387","GSM1466931","GSM1508002","GSM1466933","GSM2114337","GSM1945839","GSM1508030","GSM2114341","GSM1508007","GSM1508010","GSM1508019","GSM2029382","GSM1755407","GSM2114339","GSM1654309","GSM2029384","GSM1508008",
                "GSM1508020","GSM1654312","GSM2029385","GSM1507998","GSM2029383","GSM1755405","GSM1508028","GSM2114338","GSM1508023","GSM1654310","GSM1508021","GSM1508024","GSM1508009","GSM2029386","GSM1654311","GSM2114344","GSM1508000","GSM1508016","GSM1507997","GSM2114336","GSM2114340","GSM1508015","GSM1945838","GSM2029388","GSM2114342","GSM1508011","GSM1945837","GSM1508017","GSM1508013","GSM2114343",
                "GSM1466932","GSM1267854","GSM1294153","GSM1508022","GSM1548043","GSM1838573","GSM2033073","GSM1702299","GSM2441396","GSM2441394","GSM2441384","GSM1702300","GSM2441388","GSM2122564","GSM2122559","GSM2122562","GSM2122560","GSM2122561","GSM2122563","GSM2219677","GSM2219676","GSM2259925","GSM2259927","GSM2259926","GSM2440769","GSM2461421","GSM2461422","GSM2525732","GSM2236644","GSM2236645",
                "GSM2236646","GSM2236647","GSM2521662","GSM2521663","GSM2521664","GSM2536974","GSM2536975","GSM2536976","GSM2747880","GSM2747881","GSM2747882","GSM2747883","GSM2902360","GSM2902361","GSM2902362","GSM2844050","GSM2844051","GSM2844052","GSM2844053","GSM2151685","GSM2151686","GSM2151687","GSM2151688","GSM2151689","GSM2151690","GSM2405945","GSM2405946","GSM2405947","GSM2405948","GSM2405949",
                "GSM2405950","GSM2405951","GSM2405952","GSM2405953","GSM2405954","GSM2420187","GSM2420188","GSM2420189","GSM2464648","GSM2464649","GSM2464652","GSM2464653","GSM2464656","GSM2464657","GSM2892218","GSM2892219","GSM2892220","GSM2823388","GSM2823389","GSM2823390","GSM2823394","GSM2823395","GSM2823396","GSM2862531","GSM2862532","GSM2279160","GSM2279162")
  knockout = c("GSM2026223","GSM2026222","GSM2026232","GSM2026229","GSM2026225","GSM2026227","GSM2026224","GSM2026230","GSM2026226","GSM2026228","GSM2026221","GSM2026231","GSM1702302","GSM1702301","GSM2717435","GSM2717436","GSM2717437","GSM2717438","GSM2433705","GSM2433706","GSM2433707","GSM2433708","GSM2433709","GSM2433710","GSM2433711","GSM3100199","GSM3100200")
  
  ko = c("GSM1220020","GSM2107403","GSM1581107","GSM1580998","GSM1969203","GSM1581001","GSM1580973","GSM1581019","GSM1581011","GSM1580959","GSM1382043","GSM1581022","GSM1370702","GSM1657080","GSM1370696","GSM1581138","GSM2026223","GSM1580985","GSM1581095","GSM1580954","GSM1581102","GSM1581008","GSM1580969","GSM1581081","GSM1581025","GSM1581041","GSM1581017","GSM1581127","GSM1576148","GSM1581103","GSM1581035",
         "GSM1348967","GSM1581036","GSM1581015","GSM1580961","GSM1581111","GSM1370701","GSM1581120","GSM1665957","GSM1581010","GSM1580962","GSM1665953","GSM1581096","GSM1581075","GSM1657079","GSM1580965","GSM1581037","GSM1348965","GSM2026222","GSM1581023","GSM1580982","GSM1370704","GSM1665949","GSM1581128","GSM1657078","GSM2285686","GSM1581057","GSM1581130","GSM1581060","GSM1581026","GSM1581129",
         "GSM1382041","GSM1581134","GSM1581043","GSM1580986","GSM1581069","GSM1348966","GSM2026232","GSM1580984","GSM1581046","GSM1580976","GSM2026229","GSM1581014","GSM1581042","GSM1581137","GSM1581038","GSM1581000","GSM1581125","GSM1348954","GSM1581113","GSM2285687","GSM1665955","GSM1581005","GSM1657082","GSM1581058","GSM1581079","GSM1969206","GSM1370699","GSM1581088","GSM1665947","GSM1348964",
         "GSM1581045","GSM1969204","GSM1581100","GSM1581124","GSM1665964","GSM1581082","GSM1581059","GSM1580966","GSM1580997","GSM1580967","GSM1581009","GSM1581031","GSM1581122","GSM1581044","GSM2035725","GSM1581080","GSM1581105","GSM1581062","GSM1581104","GSM1581029","GSM1581094","GSM1581020","GSM1665960","GSM1581112","GSM1581090","GSM1581140","GSM1382042","GSM1370697","GSM1665961","GSM1581109",
         "GSM2026225","GSM1576149","GSM1581032","GSM1580963","GSM1581061","GSM1581141","GSM1580987","GSM1581108","GSM1581066","GSM1581126","GSM1581110","GSM1580968","GSM1580994","GSM1581006","GSM1581136","GSM1581083","GSM1581072","GSM1581064","GSM2026227","GSM1370695","GSM1581063","GSM1657083","GSM1370698","GSM1580978","GSM1581030","GSM1581054","GSM1581048","GSM1581073","GSM1580977","GSM1580996",
         "GSM1580991","GSM1581068","GSM1580990","GSM1581087","GSM1581040","GSM1580970","GSM1665954","GSM1581119","GSM1665962","GSM1465027","GSM1580955","GSM1665956","GSM1580981","GSM1580995","GSM1581093","GSM2026224","GSM1581047","GSM1379557","GSM1580975","GSM1580972","GSM1576144","GSM1581002","GSM1580993","GSM1370703","GSM2107402","GSM2026230","GSM2026226","GSM1581132","GSM1580964","GSM1581092",
         "GSM1581077","GSM1581028","GSM1581106","GSM2035722","GSM1581003","GSM1581049","GSM1581098","GSM1580983","GSM1580958","GSM1348955","GSM1580989","GSM1581101","GSM1581133","GSM1581004","GSM1581052","GSM2026228","GSM1581114","GSM1665951","GSM1581050","GSM1657081","GSM1580988","GSM2026221","GSM1665963","GSM1581016","GSM1348963","GSM1581065","GSM1665958","GSM1581078","GSM1581091","GSM1581121",
         "GSM1580999","GSM1581135","GSM1581089","GSM1348962","GSM1580980","GSM1580960","GSM1382044","GSM1665959","GSM1581123","GSM1581116","GSM1581033","GSM1581056","GSM1665952","GSM1581007","GSM1370700","GSM1348968","GSM1581021","GSM1581027","GSM1665948","GSM1581117","GSM1580956","GSM1581018","GSM1581074","GSM1581053","GSM1581051","GSM1465026","GSM1580979","GSM1665946","GSM1581085","GSM2285685",
         "GSM1581115","GSM1580974","GSM1576145","GSM2035723","GSM1581039","GSM1580957","GSM2035724","GSM1581086","GSM2026231","GSM1581024","GSM1665966","GSM1581099","GSM1581084","GSM1665950","GSM1581131","GSM1581012","GSM1581067","GSM1621379","GSM1939584","GSM1906571","GSM1658683","GSM1906570","GSM1906568","GSM1939589","GSM1658682","GSM1939588","GSM1658681","GSM1906569","GSM1939585","GSM1621378",
         "GSM1551899","GSM1906572","GSM1621377","GSM1906567","GSM1379556","GSM1580992","GSM1581076","GSM1581118","GSM1581071","GSM1581070","GSM1581139","GSM1581013","GSM1581097","GSM1581055","GSM1580971","GSM1581034","GSM1665965","GSM1702298","GSM1702302","GSM1702301","GSM1702297","GSM1574313","GSM1574316","GSM1574317","GSM1574314","GSM1574315","GSM2344402","GSM2344403","GSM2089695","GSM2089696",
         "GSM2089694","GSM2089693","GSM2144413","GSM2459303","GSM2493895","GSM2493897","GSM2493894","GSM2493889","GSM2493903","GSM2493892","GSM2493890","GSM2493902","GSM2493896","GSM2493901","GSM2493891","GSM2493893","GSM2543172","GSM2543186","GSM2543183","GSM2543185","GSM2543174","GSM2543176","GSM2543182","GSM2543184","GSM2543173","GSM2543175","GSM2098539","GSM2098540","GSM2098541","GSM2717435",
         "GSM2717436","GSM2717437","GSM2717438","GSM2474245","GSM2474246","GSM2474247","GSM2474248","GSM2500289","GSM2500290","GSM2500291","GSM2500292","GSM2500294","GSM2500295","GSM2500296","GSM2500297","GSM2500298","GSM2500300","GSM2500301","GSM2500303","GSM2500304","GSM2500306","GSM2500307","GSM2500308","GSM2500309","GSM2500310","GSM2515720","GSM2515721","GSM2515722","GSM2515723","GSM2864660",
         "GSM2864661","GSM2864662","GSM2864663","GSM2864664","GSM2864665","GSM2414792","GSM2104442","GSM2104443","GSM2330093","GSM2330094","GSM2330098","GSM2330100","GSM2433705","GSM2433706","GSM2433707","GSM2433708","GSM2433709","GSM2433710","GSM2433711","GSM2509513","GSM2509514","GSM2509515","GSM2509516","GSM2509517","GSM2509518","GSM2711785","GSM2711786","GSM2711787","GSM2711788","GSM2711789",
         "GSM2711790","GSM2711791","GSM2711792","GSM2711793","GSM2711794","GSM2711795","GSM2711796","GSM2786542","GSM2794657","GSM2794658","GSM2794659","GSM2794660","GSM2794661","GSM2794662","GSM2597672","GSM2597673","GSM2704295","GSM2704296","GSM2704297","GSM3021854","GSM3021855","GSM3021856","GSM3021857","GSM3021858","GSM3021859","GSM3021860","GSM3021861","GSM3118989","GSM3118990","GSM3118991",
         "GSM3118994","GSM3118995","GSM2024991","GSM2024992","GSM2024993","GSM2024994","GSM2024995","GSM2024996","GSM2024997","GSM2024998","GSM2631741","GSM2631742","GSM2631743","GSM2631744","GSM2631745","GSM2631746","GSM2631747","GSM2631748","GSM2684922","GSM2684923","GSM2684924","GSM2787815","GSM2787817","GSM2787818","GSM2787819","GSM2787820","GSM2787822","GSM2787825","GSM2787826","GSM2787827",
         "GSM2787829","GSM3100199","GSM3100200","GSM2393187","GSM2393190","GSM2393191","GSM2393195","GSM3130652","GSM3130653","GSM3130654","GSM3130655","GSM3130656","GSM3130657","GSM3130658","GSM3130659")
  
  overexpression =  c("GSM1267857","GSM1969184","GSM1548063","GSM1548068","GSM1548066","GSM1548062","GSM1378373","GSM1548069","GSM1548052","GSM1548056","GSM1548058","GSM1548065","GSM1548049","GSM1548061","GSM1548059","GSM1548051","GSM1969185","GSM1548053","GSM1548048","GSM1969186","GSM1267863","GSM1267862","GSM1548054","GSM1548050","GSM1548060","GSM1548047","GSM1548055","GSM1548057","GSM1267850","GSM1940175","GSM1969183",
                      "GSM1267851","GSM1548067","GSM1940174","GSM1381991","GSM1382000","GSM1382018","GSM1382032","GSM1382015","GSM1381984","GSM1381992","GSM1382004","GSM1382025","GSM1382017","GSM1382013","GSM1382034","GSM1382007","GSM1382020","GSM1381995","GSM1382005","GSM1381993","GSM1382024","GSM1382006","GSM1382016","GSM1382012","GSM1382010","GSM1382002","GSM1381999","GSM1382011","GSM1382014","GSM1381990",
                      "GSM1382009","GSM1382036","GSM1381989","GSM1382008","GSM1382035","GSM1382028","GSM1382023","GSM1382033","GSM1382021","GSM1381997","GSM1381987","GSM1381996","GSM1381994","GSM1382019","GSM1382001","GSM1381988","GSM1382027","GSM1381998","GSM1382022","GSM1382026","GSM1382003","GSM1548064","GSM2221672","GSM2221671","GSM2277234","GSM2277229","GSM2277232","GSM2277230","GSM2325077","GSM2325078",
                      "GSM2325079","GSM2325080","GSM2325081","GSM2325082","GSM2325083","GSM2325084","GSM2325085","GSM2325086","GSM2325089","GSM1693009","GSM1693010","GSM1693011","GSM1693012","GSM1693013","GSM1693014","GSM2627213","GSM2627214","GSM2627215","GSM2627216","GSM2627217","GSM2627218","GSM2308770","GSM2308771","GSM2944158","GSM2944160","GSM2944161","GSM2944162","GSM2944163","GSM2944164","GSM2944165",
                      "GSM2944166","GSM2944167","GSM2944168","GSM2944169")
  over_expression = c("GSM1267857","GSM1969184","GSM1548063","GSM1548068","GSM1548066","GSM1548062","GSM1378373","GSM1548069","GSM1548052","GSM1548056","GSM1548058","GSM1548065","GSM1548049","GSM1548061","GSM1548059","GSM1548051","GSM1969185","GSM1548053","GSM1548048","GSM1969186","GSM1267863","GSM1267862","GSM1548054","GSM1548050","GSM1548060","GSM1548047","GSM1548055","GSM1548057","GSM1267850","GSM1940175","GSM1969183",
                      "GSM1267851","GSM1548067","GSM1940174","GSM1381991","GSM1382000","GSM1382018","GSM1382032","GSM1382015","GSM1381984","GSM1381992","GSM1382004","GSM1382025","GSM1382017","GSM1382013","GSM1382034","GSM1382007","GSM1382020","GSM1381995","GSM1382005","GSM1381993","GSM1382024","GSM1382006","GSM1382016","GSM1382012","GSM1382010","GSM1382002","GSM1381999","GSM1382011","GSM1382014","GSM1381990",
                      "GSM1382009","GSM1382036","GSM1381989","GSM1382008","GSM1382035","GSM1382028","GSM1382023","GSM1382033","GSM1382021","GSM1381997","GSM1381987","GSM1381996","GSM1381994","GSM1382019","GSM1382001","GSM1381988","GSM1382027","GSM1381998","GSM1382022","GSM1382026","GSM1382003","GSM1548064","GSM2221672","GSM2221671","GSM2277234","GSM2277229","GSM2277232","GSM2277230","GSM2325077","GSM2325078",
                      "GSM2325079","GSM2325080","GSM2325081","GSM2325082","GSM2325083","GSM2325084","GSM2325085","GSM2325086","GSM2325089","GSM1693009","GSM1693010","GSM1693011","GSM1693012","GSM1693013","GSM1693014","GSM2627213","GSM2627214","GSM2627215","GSM2627216","GSM2627217","GSM2627218","GSM2308770","GSM2308771","GSM2944158","GSM2944160","GSM2944161","GSM2944162","GSM2944163","GSM2944164","GSM2944165",
                      "GSM2944166","GSM2944167","GSM2944168","GSM2944169")
  
  # take difference between all samples
  all_samples = unique(c(cas9, crispr, induction, knockdown, knockout, ko, overexpression, over_expression))
  
  # Retrieve information from compressed data
  samples = h5read(destination_file, "meta/Sample_geo_accession")
  tissue = h5read(destination_file, "meta/Sample_source_name_ch1")
  genes = h5read(destination_file, "meta/genes")
  # name = h5read(destination_file, 'meta/Sample_source_name_ch1') -- this is the same as tissue
  organism = h5read(destination_file, 'meta/Sample_organism_ch1')
  characteristics = h5read(destination_file, 'meta/Sample_characteristics_ch1')
  #treatment_protocol = h5read(destination_file, 'meta/Sample_treatment_protocol_ch1') -- included in geo data, but not in h5
  description = h5read(destination_file, 'meta/Sample_description')
  title = h5read(destination_file, 'meta/Sample_title')
  
  # create data.frame of all samples with associated information
  h5_metadata_df <- data.frame(organism, samples, tissue, characteristics, description, title)
  h5_metadata_df_subset <- subset(h5_metadata_df, h5_metadata_df$samples %in% all_samples)
  
  nrow(h5_metadata_df_subset)
  
  # close h5
  H5close()
  
  # remove rows with string cancer
  h5_metadata_df_subset <- h5_metadata_df_subset[!grepl('(cancer)', h5_metadata_df_subset$tissue, ignore.case = TRUE),]
  h5_metadata_df_subset <- h5_metadata_df_subset[!grepl('(cancer)', h5_metadata_df_subset$characteristics, ignore.case = TRUE),]
  h5_metadata_df_subset <- h5_metadata_df_subset[!grepl('(cancer)', h5_metadata_df_subset$description, ignore.case = TRUE),]
  
  # remove rows with string carcin
  h5_metadata_df_subset <- h5_metadata_df_subset[!grepl('(carcin)', h5_metadata_df_subset$tissue, ignore.case = TRUE),]
  h5_metadata_df_subset <- h5_metadata_df_subset[!grepl('(carcin)', h5_metadata_df_subset$characteristics, ignore.case = TRUE),]
  h5_metadata_df_subset <- h5_metadata_df_subset[!grepl('(carcin)', h5_metadata_df_subset$description, ignore.case = TRUE),]
  
  # remove rows with leukemia
  h5_metadata_df_subset <- h5_metadata_df_subset[!grepl('(leukemia)', h5_metadata_df_subset$tissue, ignore.case = TRUE),]
  h5_metadata_df_subset <- h5_metadata_df_subset[!grepl('(leukemia)', h5_metadata_df_subset$characteristics, ignore.case = TRUE),]
  h5_metadata_df_subset <- h5_metadata_df_subset[!grepl('(leukemia)', h5_metadata_df_subset$description, ignore.case = TRUE),]
  
  # remove rows with infection in characteristics -- ask
  
  # remove eomes positive nk cells in characteristics (another total rna)
  h5_metadata_df_subset <- h5_metadata_df_subset[!grepl('(eomes)', h5_metadata_df_subset$characteristics, ignore.case = TRUE),]
  
  # remove rows relating to PTSD, pre-deployment, post-deployment (ptsd studies on soldiers, rna seq)
  h5_metadata_df_subset <- h5_metadata_df_subset[!grepl('(ptsd)', h5_metadata_df_subset$characteristics, ignore.case = TRUE),]
  h5_metadata_df_subset <- h5_metadata_df_subset[!grepl('(pre-deployment)', h5_metadata_df_subset$characteristics, ignore.case = TRUE),]
  h5_metadata_df_subset <- h5_metadata_df_subset[!grepl('(post-deployment)', h5_metadata_df_subset$characteristics, ignore.case = TRUE),]
  
  # remove rows relating to neck of femure fracture hip cartilage
  h5_metadata_df_subset <- h5_metadata_df_subset[!grepl('(femure fracture hip cartilage)', h5_metadata_df_subset$tissue, ignore.case = TRUE),]
  h5_metadata_df_subset <- h5_metadata_df_subset[!grepl('(femure fracture hip cartilage)', h5_metadata_df_subset$characteristics, ignore.case = TRUE),]
  
  nrow(h5_metadata_df_subset)
  
  # de-duplicate 'tissue' column (only one replicate from each expirement will remain).
  h5_metadata_df_subset <- arrange(h5_metadata_df_subset, tissue)
  h5_metadata_df_subset <- h5_metadata_df_subset[!duplicated(h5_metadata_df_subset[,"tissue"]),]
  h5_metadata_df_subset <- h5_metadata_df_subset[!duplicated(h5_metadata_df_subset[,"characteristics"]),]
  
  print(nrow(h5_metadata_df_subset))
  
  # create vector gse with series corresponding to samples
  gse <- character(length = length(h5_metadata_df_subset$samples))
  j <- 1
  for(i in h5_metadata_df_subset$samples){
    gse[j] <- Meta(getGEO(i))$series_id
    j <- j+1
  }
  
  # add gse column to dataframe 
  human_gse <- h5_metadata_df_subset 
  human_gse[,'GSE'] <- gse
  
  # deduplicate on GSE
  human_gse <- human_gse[!duplicated(human_gse[,'GSE']),]
  
  # count rows
  print(nrow(human_gse))
  
  # get overall experimenet design from series softs
  overall_design <- character(length = length(human_gse$GSE))
  j <- 1
  for(i in human_gse$GSE){
    x <- Meta(getGEO(i, GSEMatrix = FALSE))$overall_design
    overall_design[j] <- x
    j <- j+1
  }
  
  gse_design <- human_gse['GSE']
  gse_design[, 'overall_design'] <- overall_design
  
  # write design table
  #write.csv(gse_design, 'poss_tx_datasets_seriesDesign.csv')
  #print(paste0('series design created at', getwd(), '/', 'poss_tx_datasets_seriesDesign.csv'))
  
  # write table
  #write.csv(h5_metadata_df_subset, 'human_possible_txfactor_datasets.csv')
  #print(paste0("Expression file was created at ", getwd(), "/", extracted_expression_file))
}

# extract description from GEO gse entries

# see this page re: 404 error
# https://www.biostars.org/p/377082/
