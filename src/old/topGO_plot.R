library(RSvgDevice)

setMethod("printGraph",
          signature(object = "topGOdata", result = "topGOresult",
                    firstSigNodes = "numeric", refResult = "missing"),
          function(object, result, firstSigNodes = 10, fn.prefix = "",
                   useInfo = "def", pdfSW = FALSE) {

            out.fileName <- paste(fn.prefix,  firstSigNodes, useInfo, sep = '_')              
            devSVG(file = paste(out.fileName, 'svg', sep = '.'))
            
            ## plot the graph to the specified device
            par(mai = rep(0, 4))
            gT <- showSigOfNodes(object, score(result), firstSigNodes = firstSigNodes,
                                 swPlot = FALSE, useInfo = useInfo, plotFunction = GOplot)
            plot(gT$complete.dag)
            dev.off()

            ##if(!pdfSW && .Platform$OS.type == "unix")
            ##  .ps2eps(out.fileName)
            
            cat(out.fileName, ' --- no of nodes: ', numNodes(gT$dag), '\n') 
          })



