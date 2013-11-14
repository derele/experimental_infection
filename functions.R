## selbstgeschriebene Funktion um fasta zu lesen:
read.ESeq <- function (file){
  content <- readLines(file)
  header.lines <- grep( "^>", content)
  start.lines <- header.lines+1
  end.lines <- c(header.lines[-1]-1, length(content))
  sq <- sapply(1:length(start.lines), function (i) {
    list(content[start.lines[i]:end.lines[i]])})
  names(sq) <- substr(content[header.lines],2, nchar(content[header.lines]))
  sq <- unlist(lapply(sq, paste, collapse=""))
  sq <- sq[nchar(sq)>0]
  class(sq) <- "ESeq"
  return(sq)
}

write.ESeq <- function (ESeq.obj,
                        path=paste(getwd(), "/", substitute(ESeq.obj), ".fasta", sep="")){
  write(paste(">",
              names(ESeq.obj), "\n",
              ESeq.obj, sep=""),
        file=path)
}

print.ESeq <- function(x, n=10){
  cat("* ESeq object, easy biological sequence  **\n")
  pfun <- function (x) { if(nchar(x)<10) {
    cat (x, "\n")}
  else{cat(substr(x, 1, 4), "...", "further", nchar(x)-8, "letters ...",
           substr(x, nchar(x)-3, nchar(x)),"\n")}}
  if(length(x)<n){
    sapply (x, pfun)}
  else {
    sapply(x[1:4], pfun)
    cat("...\nfurther", length(x)-8, "sequences\n...\n")
    sapply(x[(length(x)-4):length(x)], pfun)}
}

get.gc <- function (x) (nchar(gsub("A|T|N|a|t|n", "",  x))/
                        nchar(gsub("N|n", "", (x)))*100)

strReverse <- function(y) sapply(lapply(strsplit(y, NULL), rev), paste,
                                 collapse="")

revcom <- function(x){
  r <- strReverse(x)
  return(chartr("acgtryswkmbdhvnxACGTRYSWKMBDHVNX", "tgcayrswmkvhdbnxTGCAYRSWMKVHDBNX", r))
}


TOGO.all.onto <- function (o, gN, gl, g2g) {
  g <- factor(as.integer(gN %in% gl))
  names(g) <- gN
  toGO <-  new("topGOdata", ontology = o,
               allGenes = g, annot = annFUN.gene2GO,
               gene2GO = g2g)
  resultFis <- runTest(toGO, algorithm = "classic", statistic = "fisher")
  list(toGO, resultFis) ## returns a list first data then result
}
