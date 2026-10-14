library(phytools)

# read the dated tree and drop the neo-Y copy of A. juncea
grafted <- read.tree("../../divergence_dating/full_chain1/FigTree.tre")$`UTREE1=`
grafted <- drop.tip(grafted, "Aeshna_juncea_neoY")

# ladderize, then rotate nodes to this top-to-bottom species order
spp <- c("Archilestes_grandis",
         "Platycnemis_pennipes",
         "Ischnura_elegans",
         "Cordulegaster_boltonii",
         "Somatochlora_arctica",
         "Sympetrum_striolatum",
         "Orthetrum_coerulescens",
         "Pantala_flavescens",
         "Anax_junius",
         "Aeshna_juncea_neoX"
)
grafted <- ladderize(grafted)
grafted <- rotateConstr(grafted, rev(spp))

# plot leftwards with a time axis, and save the rotated tree
library(devEMF) # svg doesn't work in powerpoint
emf("figS1_tree.emf", width = 40, height = 22.05)
plot(grafted, label.offset = 10, direction = "leftwards")
axisPhylo(1)
dev.off()
write.tree(grafted, "1.odonate.synteny.nwk")

# age of the (Anax, Aeshna) node
node <- getMRCA(grafted, c("Aeshna_juncea_neoX", "Anax_junius"))
dates <- node.depth.edgelength(grafted)
age <- (max(dates) - dates[node])
age
