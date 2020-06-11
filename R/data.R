#' colorectal TCGA gene expression data with subtype annotation
#' @description The Cancer Genome Atlas level 3 RSEM genes data was retrieved
#' from Broad GDAC Firehose, doi:10.7908/C11G0KM9
#' (\url{http://gdac.broadinstitute.org/}). Samples were annotated with
#' Consensus Molecular Subtypes (Guinney 2015), Sadanandam subtypes
#' (Isella 2015) and ABSOLUTE estimated tumor purity (Aran 2015). Features were
#' annotated using \code{\link[CMScaller]{fromTo}}. Rownames are Entrez ids.
#' Included are 92 HiSeq test set samples with genes with maximum count
#' exceeding 25.
#' @references Aran D, Sirota M, Butte AJ. Systematic pan-cancer analysis of tumour purity. Nat Commun. 2015;6:8971.
#' @references Guinney J, Dienstmann R, Wang X, de Reynies A, Schlicker A, Soneson C, et al. The consensus molecular subtypes of colorectal cancer. Nat Med. 2015;21:1350-6.
#' @references Isella C, Terrasi A, Bellomo SE, Petti C, Galatola G, Muratore A, et al. Stromal contribution to the colorectal cancer transcriptome. Nat Genet [Internet]. 2015 [cited 2015 Mar 2];advance online publication. Available from: \url{http://www.nature.com/ng/journal/vaop/ncurrent/full/ng.3224.html}
#' @references TCGA. Comprehensive molecular characterization of human colon and rectal cancer. Nature. 2012;487:330-7.
#' @examples
#' library(Biobase)
#' dim(crcTCGAsubset)
"crcTCGAsubset"

#' consensus molecular subtype (CMS) templates
#' @details Colorectal cancer Consensus Molecular Subtypes (CMS) prediction
#' templates for \code{\link[CMScaller]{ntp}}. Marker genes were
#' identified using TCGA RNA-sequencing data. \code{templates$probe} refers to
#' Entrez ids.
#' @references Guinney J, Dienstmann R, Wang X, de Reynies A, Schlicker A, Soneson C, et al. The consensus molecular subtypes of colorectal cancer. Nat Med [Internet]. 2015 [cited 2015 Nov 5];advance online publication. Available from: \url{http://www.nature.com/nm/journal/vaop/ncurrent/full/nm.3967.html}
#' @references Hoshida Y. Nearest Template Prediction: A Single-Sample-Based Flexible Class Prediction with Confidence Assessment. PLoS ONE. 2010;5:e15543.
#' @examples
#' head(templates.CMS)
#' table(templates.CMS$class)
"templates.CMS.old"

#' development consensus molecular subtype (CMS) templates
#' @details Colorectal cancer Consensus Molecular Subtypes (CMS) prediction
#' templates for \code{\link[CMScaller]{ntp}}. Marker genes were
#' identified using TCGA RNA-sequencing data. \code{templates$probe} refers to
#' Entrez ids. Compared to \code{\link{templates.CMS}} this was prepared using
#' additional filters steps to further clean genes also expressed by
#' non-carcinoma cell types. Predictions are usually highly concordant with
#' \code{\link{templates.CMS}}. Under development.
#' @references Guinney J, Dienstmann R, Wang X, de Reynies A, Schlicker A, Soneson C, et al. The consensus molecular subtypes of colorectal cancer. Nat Med [Internet]. 2015 [cited 2015 Nov 5];advance online publication. Available from: \url{http://www.nature.com/nm/journal/vaop/ncurrent/full/nm.3967.html}
#' @references Hoshida Y. Nearest Template Prediction: A Single-Sample-Based Flexible Class Prediction with Confidence Assessment. PLoS ONE. 2010;5:e15543.
#' @examples
#' length(intersect(templates.CMS$symbol, templates.CMS$symbol))
#' length(setdiff(templates.CMS$symbol, templates.CMS$symbol))
"templates.CMS"


#' CRC intrinsic subtypes (CRIS) templates
#' @details Colorectal cancer CRC intrinsic subtypes (CRIS) prediction
#' templates for \code{\link[CMScaller]{ntp}} retrieved from Supplementary Table
#' 8 (Isella 2017). Suspicious (Mar-03) or legacy gene symbols were updated.
#' @references Isella C, Brundu F, Bellomo SE, Galimi F, Zanella E, Porporato R, et al. Selective analysis of cancer-cell intrinsic transcriptional traits defines novel clinically relevant subtypes of colorectal cancer. Nat Commun [Internet]. 2017 [cited 2017 Sep 11];8:ncomms15107. Available from: \url{https://www.nature.com/articles/ncomms15107}
#' @examples
#' dim(templates.CRIS)
"templates.CRIS"

#' micro-satellite instability templates
#' @details micro-satellite instability (MSI) prediction templates for
#' \code{\link[CMScaller]{ntp}}. Marker genes were identified using TCGA
#' RNA-sequencing data. \code{templates$probe} refers to Entrez ids.
#' @examples
#' head(templates.MSI)
#' table(templates.MSI$class)
"templates.MSI"

#' annotation for translating between gene identifiers
#' @description Table for translating between Ensembl, NCBI Entrez and HGNC symbol
#' human gene identifiers. Gene biotype, exonic gene length and proportion
#' GC-bases are also included. For more exhaustive annotation including gene
#' names, GO terms and so on, please refer to R Bioconductor packages
#' \code{AnnotationDbi} and \code{org.Hs.eg.db}.
#' @details Data is from GENCODE v26/GRCh38.p10. with NCBI Entrez identifiers matched using HGNC symbol against \code{\link[org.Hs.eg.db]{org.Hs.eg.db}} package.
#' @seealso \code{\link{fromTo}}, \code{\link[org.Hs.eg.db]{org.Hs.eg.db}}
#' @references Carlson M. org.Hs.eg.db: Genome wide annotation for Human. R package version 3.7.0.
#' @references Harrow J, Frankish A, Gonzalez JM, Tapanari E, Diekhans M, Kokocinski F, et al. GENCODE: The reference human genome annotation for The ENCODE Project. Genome Res. 2012;22:1760-74.
#' @examples
#' head(anno.orgHs)
#' dim(anno.orgHs)
#' colSums(apply(anno.orgHs, 2, is.na))
#' colSums(!apply(anno.orgHs[,1:4], 2, duplicated))
#' # some extremes
#' anno.orgHs[which.max(anno.orgHs$length),]
#' anno.orgHs[which.min(anno.orgHs$length),]
#' anno.orgHs[which.max(anno.orgHs$GC),]
#' anno.orgHs[which.min(anno.orgHs$GC),]
#' attributes(anno.orgHs)$genome
"anno.orgHs"

#' gene sets for exploratory gene set analysis
#' @details Gene sets from Reactome (Aug-2019) with Entrez gene identifiers.
#' @references Croft D, Mundo AF, Haw R, Milacic M, Weiser J, Wu G, et al. The Reactome pathway knowledgebase. Nucl Acids Res. 2014;42:D472-7.
#' @examples
#' head(names(geneSets.reactome))
"geneSets.reactome"

#' gene sets relevant to Consensus Molecular Subtypes
#' @description Geneset is a named list of Entrez idds. Watanabe CRC MSS/MSI,
#' Liu CDX2, Lucas HNF4A, and Servitja HNF1A were retrieved from MutSigDB C2
#' (v7.0). Gastro-Intestinal enriched genes are from the \href{http://www.proteinatlas.org/humanproteome/gastrointestinal+tract}{Human Protein Atlas}
#' Extracellular Matrix mCRC is from Naba Additional Data 2.
#' Crypt signatures are based on Merlos-Suarez Supplementary Table 5.
#' MYB signature is based on Thorner Supplementary Table 2.
#' Stromal estimate is based on Yoshihara Supplementary Data 1.
#' CRC stem cell signatures are based on de Sousa E Melo Supplementary Table 1A.
#' TGF-beta signatures are based on reanalysis of Fessler GSE79461.
#' CTNNB1/Beta-catenin signature is based on Watanabe Supplementary Table 1.
#' WNT signatures are based on Vermeulen Supplementary Table S1.
#' Retionic acid signatures are based on Duffy Supplementary Table 2.
#' Remaining are either from MutSigDB Hallmarks (v7.0) or Reactome (Aug-2019).
#' @references Croft D, Mundo AF, Haw R, Milacic M, Weiser J, Wu G, et al. The Reactome pathway knowledgebase. Nucl Acids Res. 2014;42:D472-7.
#' @references Duffy DJ, Krstic A, Halasz M, Schwarzl T, Konietzny A, Iljin K, et al. Retinoic acid and TGF-beta signalling cooperate to overcome MYCN-induced retinoid resistance. Genome Medicine. 2017;9:15.
#' @references Fessler E, Drost J, Hooff SR van, Linnekamp JF, Wang X, Jansen M, et al. TGF-beta signaling directs serrated adenomas to the mesenchymal colorectal cancer subtype. EMBO Molecular Medicine. 2016;8:745-60.
#' @references Guinney J, Dienstmann R, Wang X, de Reynies A, Schlicker A, Soneson C, et al. The consensus molecular subtypes of colorectal cancer. Nat Med. 2015;21:1350-6.
#' @references Liberzon A, Subramanian A, Pinchback R, Thorvaldsdottir H, Tamayo P, Mesirov JP. Molecular signatures database (MSigDB) 3.0. Bioinformatics. 2011;27:1739-40.
#' @references Merlos-Suarez A, Barriga FM, Jung P, Iglesias M, Cespedes MV, Rossell D, et al. The Intestinal Stem Cell Signature Identifies Colorectal Cancer Stem Cells and Predicts Disease Relapse. Cell Stem Cell. 2011;8:511-24.
#' @references Naba A, Clauser KR, Whittaker CA, Carr SA, Tanabe KK, Hynes RO. Extracellular matrix signatures of human primary metastatic colon cancers and their metastases to liver. BMC Cancer. 2014;14:518.
#' @references de Sousa E Melo F, Colak S, Buikhuisen J, Koster J, Cameron K, de Jong JH, et al. Methylation of Cancer-Stem-Cell-Associated Wnt Target Genes Predicts Poor Prognosis in Colorectal Cancer Patients. Cell Stem Cell. 2011;9:476-85.
#' @references Thorner AR, Parker JS, Hoadley KA, Perou CM. Potential Tumor Suppressor Role for the c-Myb Oncogene in Luminal Breast Cancer. PLoS One [Internet]. 2010 [cited 2017 Jan 19];5. Available from: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2951337/
#' @references Uhlen M, Fagerberg L, Hallstrom BM, Lindskog C, Oksvold P, Mardinoglu A, et al. Tissue-based map of the human proteome. Science. 2015;347:1260419.
#' @references Vermeulen L, De Sousa E Melo F, van der Heijden M, Cameron K, de Jong JH, Borovski T, et al. Wnt activity defines colon cancer stem cells and is regulated by the microenvironment. Nature Cell Biology. 2010;12:468-76.
#' @references Watanabe K, Biesinger J, Salmans ML, Roberts BS, Arthur WT, Cleary M, et al. Integrative ChIP-seq/Microarray Analysis Identifies a CTNNB1 Target Signature Enriched in Intestinal Stem Cells and Colon Cancer. PLoS One [Internet]. 2014 [cited 2017 Jan 19];9. Available from: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3961325/
#' @references Yoshihara K, Shahmoradgoli M, Martinez E, Vegesna R, Kim H, Torres-Garcia W, et al. Inferring tumour purity and stromal and immune cell admixture from expression data. Nat Commun. 2013;4:2612.
#' @examples
#' names(geneSets.CRC)
"geneSets.CRC"

#' gene sets relevant to Consensus Molecular Subtypes
#' @description Geneset is a named list of Entrez ids and is a subset of \code{\link{geneSets.CRC}}.
#' \itemize{
#' \item{MSI (Watanabe) from \href{http://software.broadinstitute.org/gsea/msigdb}{MutSigDB} C2 (v7.0)}
#' \item{DNA repair from MutSigDB Hallmark (v7.0)}
#' \item{HNF4A (Lucas) from MutSigDB C2 (v7.0)}
#' \item{MSS (Watanabe) from MutSigDB C2 (v7.0)}
#' \item{MYC from MutSigDB Hallmark (v7.0)}
#' \item{WNT signature based on Vermeulen Supplementary Table S1}
#' \item{cell cycle (E2F targets) from MutSigDB Hallmark (v7.0)}
#' \item{differentiation (gastro-intestinal markers) from
#' \href{http://www.proteinatlas.org/humanproteome/gastrointestinal+tract}{Human Protein Atlas}.}
#' \item{glycolysis from MutSigDB Hallmark (v7.0)}
#' \item{fatty acids from \href{http://www.reactome.org/}{reactome} (accessed 20190822)}
#' \item{CDX2 (Liu) from MutSigDB C2 (v7.0)}
#' \item{LGR5 stem cells Merlos-Suarez Supplementary Table 5}
#' \item{TGF-beta signature based on Fessler GEO \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE79461}{GSE79461}}
#' \item{EMT (epithelial-mesenchymal transition)from MutSigDB Hallmark (v7.0)}
#' }
#' @references Croft D, Mundo AF, Haw R, Milacic M, Weiser J, Wu G, et al. The Reactome pathway knowledgebase. Nucl Acids Res. 2014;42:D472-7.
#' @references Fessler E, Drost J, Hooff SR van, Linnekamp JF, Wang X, Jansen M, et al. TGF-beta signaling directs serrated adenomas to the mesenchymal colorectal cancer subtype. EMBO Molecular Medicine. 2016;8:745-60.
#' @references Liu T, Zhang X, So C-K, Wang S, Wang P, Yan L, et al. Regulation of Cdx2 expression by promoter methylation, and effects of Cdx2 transfection on morphology and gene expression of human esophageal epithelial cells. Carcinogenesis. 2007;28:488-96.
#' @references Lucas B, Grigo K, Erdmann S, Lausen J, Klein-Hitpass L, Ryffel GU. HNF4alpha reduces proliferation of kidney cells and affects genes deregulated in renal cell carcinoma. Oncogene. 2005;24:6418-31.
#' @references Guinney J, Dienstmann R, Wang X, de Reynies A, Schlicker A, Soneson C, et al. The consensus molecular subtypes of colorectal cancer. Nat Med. 2015;21:1350-6.
#' @references Liberzon A, Subramanian A, Pinchback R, Thorvaldsdottir H, Tamayo P, Mesirov JP. Molecular signatures database (MSigDB) 3.0. Bioinformatics. 2011;27:1739-40.
#' @references Merlos-Suarez A, Barriga FM, Jung P, Iglesias M, Cespedes MV, Rossell D, et al. The Intestinal Stem Cell Signature Identifies Colorectal Cancer Stem Cells and Predicts Disease Relapse. Cell Stem Cell. 2011;8:511-24.
#' @references Uhlen M, Fagerberg L, Hallstrom BM, Lindskog C, Oksvold P, Mardinoglu A, et al. Tissue-based map of the human proteome. Science. 2015;347:1260419.
#' @references Vermeulen L, De Sousa E Melo F, van der Heijden M, Cameron K, de Jong JH, Borovski T, et al. Wnt activity defines colon cancer stem cells and is regulated by the microenvironment. Nature Cell Biology. 2010;12:468-76.
#' @references Watanabe T, Kobunai T, Toda E, Yamamoto Y, Kanazawa T, Kazama Y, et al. Distal colorectal cancers with microsatellite instability (MSI) display distinct gene expression profiles that are different from proximal MSI cancers. Cancer Res. 2006;66:9804-8.
#' @examples
#' names(geneSets.CMS)
"geneSets.CMS"
